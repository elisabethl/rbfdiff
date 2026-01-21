close all
clear all
setPaths;
%
% An RBF-FD example for solving the Poisson equation. Collocation and LS
% unfitted or fitted methods.
%
dim = 2;                            % dim = 1,2 or 3
display = 1;                        % Plot solution
geom = 'ball';                      % ball or cube
mode = 'unfitted';                  % fitted, unfitted or collocation
bcMode = 'weak';                    % strong or weak imposition of boundary conditions (only relevant for fitted)
scaling = 1;                        % Include scaling of the unfitted LS problem
mvCentres = 1;                      % Option to have a Y point on top of all X points inside the domain
q = 3;                              % Oversampling
N = 15;                             % Number of center points (X) in each patch
P = 100;                            % Number of patches

ep = 0.1;                             % For 'phs': order of spline, 'mq', 'gs', 'iq', 'rbfqr': shape parameter, 'w2', 'bump': radius
phi = 'rbfqr';                        % Choice of basis 'phs', 'mq', 'gs', 'iq', 'rbfqr', 'w2', 'bmp'

psi = 'w2';                % Weight function: w2 or bmp
pdeg = -1;                           % Polynomial extension, not relevant for 'rbfqr'
del = 0.2;                          % Overlap between patches
%
% Place P patches and M evaluation points in geom with centre C and radius R
%
C = zeros(1,dim);
if strcmp(geom,'cube')
    R = 1*(dim)^(1/2);
elseif strcmp(geom,'ball')
    R = 1;
else
    disp("Requested geometry not implemented yet");
    return;
end
%
% Get patches  
%
ptch = getPtch(geom,P,C,R,del);
plotPtch(ptch,geom,C,R);
P = length(ptch.R);
%
% Get center points (X)
%
if strcmp(mode,"unfitted")
    dataX.nodes = [];            % Centers
    for i = 1:P
        [ptch.xc(i)] = getPts("ball",N,0,ptch.C(i,:),ptch.R(i),"unfitted",0);
        xCglobalId{i} = [(i-1)*N+1:i*N]';
        dataX.nodes = [dataX.nodes; ptch.xc(i).nodes];
    end
    xc = dataX.nodes;
    [ptch.xc.globalId] = xCglobalId{:};
    if strcmp(geom,"ball")
        dataX.inner = find(sqrt(sum((xc-C).^2,2))<=R);
    elseif strcmp(geom,"cube")
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];  
        dataX.inner = find(all(abs(xc-C)<=dimLCoeff(dim)*R,2));
    end
    dataX.outer = setdiff(1:size(xc,1),dataX.inner);
    dataX.bnd = [];
elseif strcmp(mode,"collocation")
    dataX = getPts(geom,N*P,0,C,R,"fitted",0);
    xc = dataX.nodes;
    for i = 1:P
        ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
        ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
    end
elseif strcmp(mode,"fitted") 
    dataX = getPts(geom,N*P,0,C,R,"fitted",0);
    xc = dataX.nodes;
    for i = 1:P
        ptch.xc(i).globalId = find(sqrt(sum((xc - ptch.C(i,:)).^2,2)) <= ptch.R(i));
        ptch.xc(i).nodes = xc(ptch.xc(i).globalId,:);
    end
end
%
% Get evaluation points (Y)
% 
q = max(q*double(~strcmp(mode,"collocation")),1);
M = ceil(N*P*q);  
if ~strcmp(mode,"collocation")
    dataY = getPts(geom,M,0,C,R,"fitted",0);
    %
    % Move evaluation points inside (and on the boundary) to the closest center point.
    %
    if mvCentres && ~strcmp(mode,"collocation")
        dataY = movePts(dataX,dataY);
    end
    xe = dataY.nodes;
    for i = 1:P
        ptch.xe(i).globalId = find(sqrt(sum((xe - ptch.C(i,:)).^2,2)) <= ptch.R(i));
        ptch.xe(i).nodes = xe(ptch.xe(i).globalId,:);
    end
else
    dataY = dataX;
    xe = xc;        % Centers
    ptch.xe = ptch.xc;
end

%
% Constructing global RBF-PUM approximation to evaluation and Laplace operators M x N
%
Eglobal = spalloc(M,P*N,M*N);
Lglobal = spalloc(M,P*N,M*N);
[w] = weights(psi,1.5,ptch);
for i = 1:P
    Psi = RBFInterpMat(phi,pdeg,ep,ptch.xc(i).nodes,ptch.C(i,:),ptch.R(i));
    E = RBFDiffMat(0,Psi,ptch.xe(i).nodes);
    B = RBFDiffMat(1,Psi,ptch.xe(i).nodes);
    L = RBFDiffMat(1.5,Psi,ptch.xe(i).nodes);
    Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.f'.*E + Eglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
    for d = 1:dim
        Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) =  2.*w{i}.grad{d}'.*B{d} + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
    end
    Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId) = w{i}.L'.*E + w{i}.f'.*L + Lglobal(ptch.xe(i).globalId,ptch.xc(i).globalId);
end
L = Lglobal(dataY.inner,:);
B = Eglobal(dataY.bnd,:);
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xe(dataY.inner,1)); fun(xe(dataY.bnd,1))];
    uExact = fun(xe(:,1));
    ucExact = fun(xc(:,1));
    lapAnalytic = lapFun(xe(dataY.inner,1));
    bndAnalytic = [0; 0];
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(dataY.inner,1),xe(dataY.inner,2)); fun(xe(dataY.bnd,1),xe(dataY.bnd,2))];
    uExact = fun(xe(:,1),xe(:,2));
    ucExact = fun(xc(:,1),xc(:,2));
    lapAnalytic = lapFun(xe(dataY.inner,1),xe(dataY.inner,2));
    bndAnalytic = fun(xe(dataY.bnd,1),xe(dataY.bnd,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xe(dataY.inner,1),xe(dataY.inner,2),xe(dataY.inner,3)); fun(xe(dataY.bnd,1),xe(dataY.bnd,2),xe(dataY.bnd,3))];
    uExact = fun(xe(:,1),xe(:,2),xe(:,3));
    ucExact = fun(xc(:,1),xc(:,2),xc(:,3));
    lapAnalytic = lapFun(xe(dataY.inner,1),xe(dataY.inner,2),xe(dataY.inner,3));
    bndAnalytic = fun(xe(dataY.bnd,1),xe(dataY.bnd,2),xe(dataY.bnd,3));
end 
%
% LS problem scaling
%
if scaling && (strcmp(mode,"unfitted") || (strcmp(mode,"fitted") && strcmp(bcMode,"weak")))
    [~,dist] = knnsearch(xc,xc,'K',2);
    h = max(dist(:,2));
    Lscale = sqrt(dataY.Vol/length(dataY.inner));
    Bscale = sqrt(dataY.Area/length(dataY.bnd))*(h.^-1.5);
    L = Lscale.*L; B = Bscale.*B;
    F(dataY.inner) = Lscale.*F(dataY.inner);
    F(dataY.bnd) = Bscale.*F(dataY.bnd);
end
%
% Impose boundary conditions strongly for fitted method
%
if strcmp(mode,"fitted") && strcmp(bcMode,"strong")
    Fmod = Lglobal(:,dataX.bnd)*ucExact(dataX.bnd);
    F = F-Fmod;
    L = Lglobal(:,dataX.inner);
    B = zeros(0,N);
end
%
% Solution on centre points, evaluated on Y set
%
A = [L; B];
% Ltest = rank(full(L))
% Btest = rank(full(B))
u = A\F;
%
% Fix operators to compute error measures
%
if strcmp(mode,"fitted") && strcmp(bcMode,"strong")
    u = [u; ucExact(dataX.bnd)];
end
ue = Eglobal*u;
l2Error = norm(abs(ue-uExact),2)/norm(uExact,2);
%
% Displaying output 
%
if display
    %
    % Plotting solution and patches
    %
    plotPtch(ptch,geom,C,R)
    hold on
    if dim == 1
        evalPtPlot = plot(xe,zeros(size(xe,1),1),'b.');
        centerPtPlot = plot(xc,zeros(size(xc,1),1),'rx');
    elseif dim == 2
        evalPtPlot = plot(xe(:,1),xe(:,2),'b.');
        centerPtPlot = plot(xc(:,1),xc(:,2),'rx');
    elseif dim == 3
        evalPtPlot = plot3(xe(:,1),xe(:,2),xe(:,3),'b.');
        hold on
        centerPtPlot = plot3(xc(:,1),xc(:,2),xc(:,3),'rx');
    end
    legend([evalPtPlot centerPtPlot],'Evaluation points','Center points','FontSize',18,'Interpreter','latex');
    
    plotSolution(ue,uExact,xe,dataY.bnd,dim,geom);
    %
    % Operator and solution l2 errors
    %
    lapNumeric = Lglobal(dataY.inner,:)*ucExact;
    evalNumeric = Eglobal(dataY.bnd,:)*ucExact;
    laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
    bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
    disp(['PDE error = ', num2str(l2Error)]);
    disp(['Boundary Op error = ', num2str(bndError)]);
    disp(['Laplace Op error = ', num2str(laplaceError)]);
end
%
% Solution and error plotting routines
%
function [] = plotSolution(uNumeric,uAnalytic,x,idB,dim,geom)
    error = uNumeric-uAnalytic;
    if dim == 1
        [x,id] = sort(x);
        figure()
        plot(x,uAnalytic(id),'r-'); 
        hold on;
        plot(x,uNumeric(id),'bo');
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        legend('$u_E$','$u_N$','FontSize',18,'Interpreter','latex');
        grid on
    
        figure()
        plot(x,error(id),'b-.');
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("$e_{rel}$","Interpreter","latex","FontSize",24,'Rotation',0)
        grid on
    elseif dim == 2
        T = delaunay(x(:,1),x(:,2));

        figure()
        G=trisurf(T,x(:,1),x(:,2),uAnalytic);
        hold on
        if strcmp(geom,"ball")
            plot(x(idB,1),x(idB,2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
        set(G,'EdgeColor','none')
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),uNumeric);
        hold on
        if strcmp(geom,"ball")
            plot(x(idB,1),x(idB,2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        set(G,'EdgeColor','none')
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),error);
        hold on
        if strcmp(geom,"ball")
            plot(x(idB,1),x(idB,2),'k-',"LineWidth",1.5)
        elseif strcmp(geom,"cube")
            plot([max(x(idB,1)) max(x(idB,1)) min(x(idB,1)) min(x(idB,1)) max(x(idB,1))],...
                 [max(x(idB,2)) min(x(idB,2)) min(x(idB,2)) max(x(idB,2)) max(x(idB,2))],'k-',"LineWidth",1.5)
        end
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u-u_E$$","Interpreter","latex","FontSize",24)
        set(G,'EdgeColor','none')
        shading interp
    elseif dim == 3 && ~strcmp(geom,"cube")
        
        bndPtCloud = pointCloud(x(idB,:));
        surfMesh = pc2surfacemesh(bndPtCloud,"ball-pivot");
        [~,~,idRes] = intersect(surfMesh.Vertices,x(idB,:),'stable','rows');
        x = surfMesh.Vertices;
        T = surfMesh.Faces;
        
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uAnalytic(idB(idRes)));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
        shading interp
        
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uNumeric(idB(idRes)));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
        shading interp
    
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),error(idB(idRes)));
        axis equal
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("z","Interpreter","latex","FontSize",24,"Rotation",0)
        set(G,'EdgeColor','none')
        cb = colorbar;
        ylabel(cb,"$$u-u_E$$","Interpreter","latex","FontSize",24,'Rotation',90)
        shading interp
    else
        disp("No plotting available in 3D for this geometry")
    end
end
%
% Patch plotting routine
%
function [] = plotPtch(ptch,geom,C,R)
    theta = linspace(0,2*pi,100)';
    dim = size(C,2);
    figure()
    %
    % Plot patches
    %
    if dim == 1
        altern = 0;
        for i = 1:size(ptch.C,1)
            plot(ptch.C(i,1),altern,'ro');
            hold on
            x = [ptch.C(i,1) + ptch.R(i); ptch.C(i,:) - ptch.R(i)];
            plot(x(:,1),altern.*ones(2,1),'k-')
            altern = mod(i,2)*1e-2;
        end
    elseif dim == 2
        for i = 1:size(ptch.C,1)
            plot(ptch.C(i,1),ptch.C(i,2),'ro');
            hold on
            x = ptch.C(i,:) + [ptch.R(i)*cos(theta), ptch.R(i)*sin(theta)];
            plot(x(:,1),x(:,2),'k-')
        end
    elseif dim == 3
        [x,y,z] = sphere(10);
        for i = 1:size(ptch.C,1)
            plot3(ptch.C(i,1),ptch.C(i,2),ptch.C(i,3),'ro');
            hold on
            surf(x.*ptch.R(i)+ptch.C(i,1),y.*ptch.R(i)+ptch.C(i,2),z.*ptch.R(i)+ptch.C(i,3),...
                'FaceColor','none');
        end
    end
    %
    % Plot geometry frame
    %
    if strcmp(geom,"cube")
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];   
        a = dimLCoeff(dim)*R;
        if dim == 1
            x = [C + a; ...
                 C - a];
            plot(x(:,1),zeros(size(x,1)),'ko');
        elseif dim==2
            x = [C + [a a]; ...
                 C + [a -a]; ...
                 C + [-a -a]; ...
                 C + [-a a]; ...
                 C + [a a]];
            plot(x(:,1),x(:,2),'k-');
        elseif dim==3
            x = [C + [a a a]; ...
                 C + [a -a a]; ...
                 C + [-a -a a]; ...
                 C + [-a a a]; ...
                 C + [a a a]];
            x = [x; [x(:,1), x(:,2), x(:,3)-2*a]];
            plot3(x(:,1),x(:,2),x(:,3),'k-','LineWidth',2);
            for i = 1:size(x,1)/2
                xPlot = [x(i,:); x(i+size(x,1)/2,:)];
                plot3(xPlot(:,1),xPlot(:,2),xPlot(:,3),'k-','LineWidth',2);
            end
        end
    elseif strcmp(geom,"ball")
        if dim == 1
            x = [C + R; ...
                 C - R];
            plot(x(:,1),zeros(size(x,1)),'ko');
        elseif dim == 2
            theta = linspace(0,2*pi,1000);
            plot(R.*cos(theta),R.*sin(theta),'k-')
        elseif dim == 3
            [x,y,z] = sphere(20);
            surf(x.*R+C(1),y.*R+C(2),z.*R+C(3),'FaceColor','none','LineWidth',2);
        end
    end
    axis equal
end
%
% Point generation routine
%
function data = getPts(geom,N,n,C,R,mode,extCoeff)
    dim = size(C,2);
    if strcmp(geom,'ball')
        dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];       % Ratios between rectangle/circle, cube/sphere area and volume
        dimACoeff = [1 pi (4/3)*pi];                % constant for size of domain
        dimPCoeff = [1 2*pi 4*pi];                  % constant for computing size of boundary
        h = R/(0.5*N^(1/dim) - extCoeff*n^(1/dim)); % approximate fill distance
        xB = zeros(0,dim);
        %
        % If using the unfitted method, no points on the boundary
        %
        if ~strcmp(mode,'unfitted')
            % Boundary point generation
            if dim == 1
                xB = [-R, R]';
            elseif dim == 2
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*h^(dim-1)));
                xB = [R*cos(linspace(-pi,pi-(2*pi)/(Nb+1),Nb)'), R*sin(linspace(-pi,pi-(2*pi)/(Nb+1),Nb)')];
            elseif dim == 3
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*(0.5*h)^(dim-1)));
                % Fibonacci points on sphere
                ratio = 1+sqrt(5);
                ind = [0:Nb-1]' + 0.5;
                theta = pi*ratio*ind;
                fi = acos(1-(2*ind)/(Nb));
                xB = [R.*sin(fi).*cos(theta), R.*sin(fi).*sin(theta), R.*cos(fi)];
            end
        end
        Nb = length(xB);
        % Interior point generation
        Ncube = ceil(dimRat(dim)*(N-Nb));
        x = 2*(R+n^(1/dim)*h*extCoeff)*(halton(Ncube,dim)-0.5);
        r2 = sqrt(sum(x.^2,2));
        pos = find(r2<=(R+n^(1/dim)*h*extCoeff));
        pos = pos(1:N-Nb);
        x = x(pos,:) + C;
        [~,ii] = sort(sum(x,2));
        x = x(ii,:);
        %
        % Organize outputs including labels for points inside, outside and on the boundary
        %
        data.nodes = [x; xB];
        data.inner = find(sqrt(sum((x-C).^2,2))<=R);
        data.outer = setdiff(1:size(x,1),data.inner)';
        data.bnd = [size(x,1) + 1:size(x,1) + size(xB,1)]';
        %
        % Domain volume and area measures used for the scaling of the LS problem
        %
        data.Vol = (dimACoeff(dim)*R^dim);
        data.Area = (dimPCoeff(dim)*R^(dim-1));
    elseif strcmp(geom,'cube')
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];          % constants for getting cube side length /2
        h = R/(0.5*N^(1/dim) - extCoeff*n^(1/dim));   % approximate fill distance
        xB = zeros(0,dim);
        %
        % If using the unfitted method, no points on the boundary
        %
        if ~strcmp(mode,'unfitted')
            % Boundary point generation
            if dim == 1
                xB = [-R, R]';
            elseif dim == 2
                Nb = ceil(((2*dimLCoeff(dim)*R)^(dim-1))/(h^(dim-1)));
                XYZLim = [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                xB = [linspace(XYZLim(1),XYZLim(2),Nb)', ones(Nb,1)*XYZLim(1);...
                      ones(Nb,1)*XYZLim(1), linspace(XYZLim(1),XYZLim(2),Nb)';...
                      linspace(XYZLim(1),XYZLim(2),Nb)', ones(Nb,1)*XYZLim(2);...
                      ones(Nb,1)*XYZLim(2), linspace(XYZLim(1),XYZLim(2),Nb)'];
                xB = unique(xB,'rows');
                xB = xB + C;
            elseif dim == 3
                NbF = ceil(((2*dimLCoeff(dim)*R)^(dim-1))/(h^(dim-1)));     % Number of points on each face
                NbE = ceil(2*dimLCoeff(dim)/h);                             % Number of points on each edge
                XYZLim = [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                xB = [linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1), ones(NbE,1)*XYZLim(1);...
                      ones(NbE,1)*XYZLim(1), linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1);...
                      linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(2), ones(NbE,1)*XYZLim(1);...
                      ones(NbE,1)*XYZLim(2), linspace(XYZLim(1),XYZLim(2),NbE)', ones(NbE,1)*XYZLim(1)];
                xB = [xB; [1,1,-1].*xB];
                Ry = eye(dim); Ry(1,1) = cos(pi/2); Ry(dim,dim) = cos(pi/2); Ry(dim,1) = -sin(pi/2); Ry(1,dim) = sin(pi/2); % Rotate 90 degrees along y-axis
                xB = [xB; (Ry*xB')'];
                xB = uniquetol(xB,1e-14,'ByRows',true);
                NbE = size(xB,1);
                Rx = eye(dim); Rx(dim-1,dim-1) = cos(pi/2); Rx(dim,dim) = cos(pi/2); Rx(dim-1,dim) = -sin(pi/2); Rx(dim,dim-1) = sin(pi/2); % Rotate 90 degrees
                xBF = [2*(dimLCoeff(dim)*R)*(halton(NbF,dim-1)-0.5), ones(NbF,1)*XYZLim(1)];
                xBF = [xBF; [1,1,-1].*xBF];
                xB = [xB; xBF];
                xB = [xB; (Rx*xBF')'];
                xB = [xB; (Ry*xBF')'];
            end
        end
        Nb = length(xB);
        % Interior point generation
        x = 2*(dimLCoeff(dim)*R+dimLCoeff(dim)*n^(1/dim)*h*extCoeff)*(halton(N-Nb,dim)-0.5);
        r2 = sqrt(sum(x.^2,2));
        x = x + C;
        %
        % Organize outputs including labels for points inside, outside and on the boundary
        %
        data.nodes = [x; xB];
        data.inner = find(all(abs(x-C)<=dimLCoeff(dim)*R,2));
        data.outer = setdiff(1:size(x,1),data.inner)';
        data.bnd = [size(x,1) + 1:size(x,1) + size(xB,1)]';
        %
        % Domain volume and area measures used for the scaling of the LS problem
        %
        data.Vol = (2*dimLCoeff(dim)*R)^dim;
        data.Area = 2*dim*(2*dimLCoeff(dim)*R)^(dim-1);
    end
end
%
% Patch generation routine
%
function ptch = getPtch(geom,P,C,R,del)
    dim = size(C,2);
    if strcmp(geom,"ball")
        dimRat = [1 4/pi 8/(4*pi/3)];
        P = ceil(dimRat(dim)*P);
        Rsq = R*(dim)^(1/2);
    elseif strcmp(geom,"cube")
        Rsq = R;
    end
    dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];              % constants for getting cube side length /2
    ptch.R = (2*Rsq)/((2-del)*ceil(P^(1/dim)) - del); % Ensure overlap in cartesian domain (Along diagonal)
    if dim == 1
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        ptch.C = Cx(:);
    elseif dim == 2
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy] = meshgrid(Cx,Cy);
        ptch.C = [Cx(:),Cy(:)];
    elseif dim == 3
        Cx = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cy = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        Cz = linspace(-1+dimLCoeff(dim)*(1-del)*ptch.R,1-dimLCoeff(dim)*(1-del)*ptch.R,ceil(P^(1/dim)));
        [Cx,Cy,Cz] = meshgrid(Cx,Cy,Cz);
        ptch.C = [Cx(:),Cy(:),Cz(:)];
    end
    ptch.R = ptch.R*ones(size(ptch.C,1),1);
    if strcmp(geom,"ball")
        tol = 1e-6;
        % Get closest patch neighbours based on distance
        [idNei,distNei] = rangesearch(ptch.C,ptch.C,(2*ptch.R(1)-ptch.R(1)*del) + ptch.R(1)*tol);
        % Find points of intersection between patches (2 patches in 2D and 4 in 3D)
        edgePts = [];
        for i = 1:size(idNei,1)
            neiTemp = idNei{i};
            distTemp = distNei{i};
            % All possible patch combinations in neighbourhood, groups of 2 in 2D and 4 in 3D.
            ptchComb = [nchoosek(neiTemp(2:end),2^(dim-1) - 1)];
            ptchComb = [ones(size(ptchComb,1),1).*neiTemp(1), ptchComb];
            for j = 1:size(ptchComb,1)
                midPt = sum(ptch.C(ptchComb(j,:),:),1)/size(ptchComb,2);
                tempDist = sqrt(sum((ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:)).^2,2));
                tempDir = (ptch.C(ptchComb(j,2:end),:) - ptch.C(ptchComb(j,1),:))./tempDist;
                % Computing directions along which the edge point lies (starting from midPt)
                if dim == 2
                    t = (-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:))/norm((-ptch.C(neiTemp(1),:) + ptch.C(neiTemp(j+1),:)),2);
                    t = [t(2), -t(1)];
                    t2 = -t;
                    keep = 1;
                else
                    t = cross(tempDir(1,:),tempDir(2,:));
                    keep = ~dot(t,tempDir(end,:));
                    t2 = cross(tempDir(1,:),-tempDir(2,:));
                end
                if keep % In 3D some combinations of patch centres do not form a plane
                    tempL = sqrt((ptch.R(1)^2)-(max(tempDist)/2)^2);
                    edgePts = [edgePts; tempL*t + midPt];
                    edgePts = [edgePts; tempL*t2 + midPt];                     
                end                    
            end
        end
        edgePts= uniquetol(edgePts,tol,'ByRows',true);
        % Make sure that points of intersection that are close to or inside
        % the domain are fully covered. This guarantees no part of the
        % boundary is left uncovered.
        idEdgeIn = find(sum((edgePts - C).^2,2)<=(R+del.*ptch.R(1)).^2);
        edgePtsIn = edgePts(idEdgeIn,:);
        [~,ptchOrder] = sort(sum((ptch.C - C).^2,2),'ascend');
        idBall = [];
        for i = 1:size(edgePtsIn,1)
            ptPtchList = sqrt(sum((edgePtsIn(i,:) - ptch.C(ptchOrder,:)).^2,2)) < ptch.R(ptchOrder) - ptch.R(ptchOrder).*tol;
            idCand = find(ptPtchList);
            idBall = [idBall ptchOrder(idCand(1))];
        end
        % idBall = find(sum(ptPtchList(idEdgeIn,:),1));
        idBall = unique(idBall);
        ptch.C = ptch.C(idBall,:);
        ptch.R = ptch.R(idBall,:);
    end
end
%
% Make sure not to move two evaluation points to the same center point. All
% center points should be moved to. Move interior and then boundary points
% to avoid changing boundary
%
function dataY = movePts(dataX,dataY)
    % Move interior points
    xcIn = dataX.nodes(dataX.inner,:);
    xeIn = dataY.nodes(dataY.inner,:);
    i = 1;
    idY = [1:size(xeIn,1)]'; % All interior evaluation (Y) points
    idX = [1:size(xcIn,1)]'; % All interior center (X) points
    %
    % Ensure the same Y point is not moved twice and obly consider X
    % points that have not been treated
    %
    iXdone = [];             
    iYdone = [];
    while length(iXdone) ~= size(xcIn,1) % ensure all interior center points are treated
        idYmv = knnsearch(xeIn(idY,:),xcIn(idX,:),'k',1);
        [unIdYmv,iX,~] = unique(idYmv,'first');
        xeIn(idY(unIdYmv),:) = xcIn(idX(iX),:);
        iXdone = [iXdone; idX(iX)];
        iYdone = [iYdone; idY(unIdYmv)];
        i = i + 1;
        idY = setdiff(idY,iYdone);
        idX = setdiff(idX,iXdone);
    end
    dataY.nodes(dataY.inner,:) = xeIn;
    %
    % Move boundary points for fitted method only
    %
    if ~isempty(dataX.bnd)
        xcBnd = dataX.nodes(dataX.bnd,:);
        xeBnd = dataY.nodes(dataY.bnd,:);
        i = 1;
        idY = [1:size(xeBnd,1)]'; % All interior evaluation (Y) points
        idX = [1:size(xcBnd,1)]'; % All interior center (X) points
        %
        % Ensure the same Y point is not moved twice and obly consider X
        % points that have not been treated
        %
        iXdone = [];             
        iYdone = [];
        while length(iXdone) ~= size(xcBnd,1) % ensure all interior center points are treated
            idYmv = knnsearch(xeBnd(idY,:),xcBnd(idX,:),'k',1);
            [unIdYmv,iX,~] = unique(idYmv,'first');
            xeBnd(idY(unIdYmv),:) = xcBnd(idX(iX),:);
            iXdone = [iXdone; idX(iX)];
            iYdone = [iYdone; idY(unIdYmv)];
            i = i + 1;
            idY = setdiff(idY,iYdone);
            idX = setdiff(idX,iXdone);
        end
        dataY.nodes(dataY.bnd,:) = xeBnd;
    end
end
