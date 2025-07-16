close all
clear all
setPaths;
%
% An RBF-FD example for solving the Poisson equation. Collocation and
% unfitted, fitted methods.
%
dim = 2;                            % dim = 1,2 or 3
display = 1;                        % Plot solution
geom = 'cube';                      % ball or cube
mode = 'unfitted';                  % fitted, unfitted or collocation
scaling = 1;                        % Include scaling of the unfitted LS problem
mvCentres = 0;                      % Option to have a Y point on top of all X points inside the domain
q = 2;                              % Oversampling
N = 500;                            % Number of center points
ep = 1;                             % Not relevant for 'r3' basis
phi = 'r3';                         % Choice of basis 'r3', 'mq', 'gs', 'iq', 'rbfqr'
pdeg = 2;                           % Polynomial extension, not relevant for 'rbfqr'
rbfDeg = 4;                         % Relevant when there is no polynomial extension
extCoeff = 0.5;                     % Extension size (in % of stencil), relevant for unfitted method
if pdeg == -1
    n = nchoosek(rbfDeg+dim-1,dim); % Order in smooth RBF case
else
    n = 2*nchoosek(pdeg+dim,dim);   % Stencil size
end
extCoeff = extCoeff*(strcmp(mode,"unfitted"));
%
% Place N centre points and M evaluation points in geom with centre C and radius R
%
C = zeros(1,dim);
R = sqrt(2);%1;
[xc,xcB,V,Area] = getPts(geom,N,n,C,R,mode,extCoeff);
xcAll = [xc; xcB];                                      % Interior and boundary points 
idXin = 1:size(xc,1);                                   % Interior pts index
idXbnd = size(xc,1)+1:size(xcAll,1);                    % Boundary pts index
% 
% For collocation, evaluation and centre points match
%
if strcmp(mode,"collocation")
    M = N;
    xe = xc; xeB = xcB;
else
    M = q*N;
    [xe,xeB] = getPts(geom,M,n,C,R,"fitted",0);
    %
    % Move evaluation points inside to the closest center point
    %
    if mvCentres
        xcIn = xc(sqrt(sum(xc.^2,2))<=R,:);
        xeTemp = xe;
        for i = 1:length(xcIn)
            idMove = knnsearch(xeTemp,xcIn(i,:),'K',1);
            xe(idMove,:) = xcIn(i,:);
            xeTemp(idMove,:) = 2*R.*ones(1,dim);
        end
        idMove = knnsearch(xe,xcIn,'K',1);
        xe(idMove,:) = xcIn;
    end
end
xeAll = [xe; xeB];                      % Interior and boundary points 
idYin = 1:size(xe,1);                   % Interior pts index
idYbnd = size(xe,1)+1:size(xeAll,1);    % Boundary pts index
%
% Ensure center points are not too far from boundary and make evaluation point - stencil list 
%
if strcmp(mode,"unfitted")
    xcAll = xcAll(unique(knnsearch(xcAll,xeAll,'K',ceil(extCoeff*n))),:);
end
N = length(xcAll);
ptStencilList = knnsearch(xcAll,xeAll,'K',1);
relXc = unique(ptStencilList); % stencil centres that are used directly
%
% Constructing global LS-RBF-FD approximation to evaluation and Laplace operators M x N
%
Eglobal = spalloc(M,N,M*n);
Lglobal = spalloc(M,N,M*n);
for i = 1:length(relXc)
    [idX,~] = knnsearch(xcAll,xcAll(relXc(i),:),'K',n);
    xcLoc = xcAll(idX,:); % stencil
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    idY = find(ptStencilList==relXc(i));
    E = RBFDiffMat(0,Psi,xeAll(idY,:));
    L = RBFDiffMat(1.5,Psi,xeAll(idY,:));
    
    Eglobal(idY,idX) = E + Eglobal(idY,idX);
    Lglobal(idY,idX) = L + Lglobal(idY,idX);
end
L = Lglobal(idYin,:);
B = Eglobal(idYbnd,:);
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xe(:,1)); fun(xeB(:,1))];
    uExact = fun(xeAll(:,1));
    ucExact = fun(xcAll(:,1));
    lapAnalytic = lapFun(xe(:,1));
    bndAnalytic = [0; 0];
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(:,1),xe(:,2)); fun(xeB(:,1),xeB(:,2))];
    uExact = fun(xeAll(:,1),xeAll(:,2));
    ucExact = fun(xcAll(:,1),xcAll(:,2));
    lapAnalytic = lapFun(xe(:,1),xe(:,2));
    bndAnalytic = fun(xeB(:,1),xeB(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xe(:,1),xe(:,2),xe(:,3)); fun(xeB(:,1),xeB(:,2),xeB(:,3))];
    uExact = fun(xeAll(:,1),xeAll(:,2),xeAll(:,3));
    ucExact = fun(xcAll(:,1),xcAll(:,2),xcAll(:,3));
    lapAnalytic = lapFun(xe(:,1),xe(:,2),xe(:,3));
    bndAnalytic = fun(xeB(:,1),xeB(:,2),xeB(:,3));
end 
%
% LS problem scaling
%
if scaling && strcmp(mode,"unfitted")
    [~,dist] = knnsearch(xcAll,xcAll,'K',2);
    h = max(dist(:,2));
    Lscale = sqrt(V/length(idYin));
    Bscale = sqrt(Area/length(idYbnd))*(h.^-1.5);
    L = Lscale.*L; B = Bscale.*B;
    F(idYin) = Lscale.*F(idYin);
    F(idYbnd) = Bscale.*F(idYbnd);
end

if strcmp(mode,"fitted")
    Fmod = Lglobal(:,idXbnd)*ucExact(idXbnd);
    F = F-Fmod;
    L = Lglobal(:,idXin);
    B = zeros(0,N);
end
%
% Solution on centre points, evaluated on Y set
%
A = [L; B];
u = A\F;
%
% Fix operators to compute error measures
%
if strcmp(mode,"fitted")
    u = [u; ucExact(idXbnd)];
end
ue = Eglobal*u;
%
% Displaying output 
%
if display
    % Plotting 
    plotSolution(ue,uExact,xeAll,xeB,dim);
    %
    % Operator and solution l2 errors
    %
    lapNumeric = Lglobal(idYin,:)*ucExact;
    evalNumeric = Eglobal(idYbnd,:)*ucExact;
    l2Error = norm(abs(ue-uExact),2)/norm(uExact,2);
    laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
    bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
    disp(['PDE error = ', num2str(l2Error)]);
    disp(['Boundary Op error = ', num2str(bndError)]);
    disp(['Laplace Op error = ', num2str(laplaceError)]);
end
%
% Plotting routines
%
function [] = plotSolution(uNumeric,uAnalytic,x,xB,dim)
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
        plot(xB(:,1),xB(:,2),'k-',"LineWidth",1.5)
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
        plot(xB(:,1),xB(:,2),'k-',"LineWidth",1.5)
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
        plot(xB(:,1),xB(:,2),'k-',"LineWidth",1.5)
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("y","Interpreter","latex","FontSize",24)
        zlabel("$$u-u_E$$","Interpreter","latex","FontSize",24)
        set(G,'EdgeColor','none')
        shading interp
    elseif dim == 3

        bndPtCloud = pointCloud(xB);
        surfMesh = pc2surfacemesh(bndPtCloud,"ball-pivot");
        x = surfMesh.Vertices;
        T = surfMesh.Faces;
        
        figure()
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uAnalytic(end-length(xB)+1:end));
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
        G=trisurf(T,x(:,1),x(:,2),x(:,3),uNumeric(end-length(xB)+1:end));
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
        G=trisurf(T,x(:,1),x(:,2),x(:,3),error(end-length(xB)+1:end));
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
    end
end
%
% Point generation routine
%
function [x,xB,Vol,Area] = getPts(geom,N,n,C,R,mode,extCoeff)
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
        %
        % Domain volume and area measures used for the scaling of the LS
        % problem
        %
        Vol = (dimACoeff(dim)*R^dim);
        Area = (dimPCoeff(dim)*R^(dim-1));
    end
    
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
        %
        % Domain volume and area measures used for the scaling of the LS
        % problem
        %
        Vol = (dimACoeff(dim)*R^dim);
        Area = (dimPCoeff(dim)*R^(dim-1));
    elseif strcmp(geom,'cube')
        dimACoeff = [1 2 8*sqrt(3)/9];                % constant for size of domain
        dimPCoeff = [1 4*sqrt(2) 8];                  % constant for computing size of boundary
        dimLCoeff = [1 sqrt(2)/2 sqrt(3)/3];          % constant for computing size of boundary
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
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*h^(dim-1)));
                xLim = C(1) + [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                yLim = C(2) + [-dimLCoeff(dim)*R dimLCoeff(dim)*R];
                xB = [linspace(xLim(1),xLim(2),floor(Nb/4))', ones(floor(Nb/4),1)*yLim(1);...
                      ones(floor(Nb/4),1)*xLim(1), linspace(yLim(1),yLim(2),floor(Nb/4))';...
                      linspace(xLim(1),xLim(2),floor(Nb/4))', ones(floor(Nb/4),1)*yLim(2);...
                      ones(floor(Nb/4),1)*xLim(2), linspace(yLim(1),yLim(2),floor(Nb/4))'];
                xB = unique(xB,'rows');
            elseif dim == 3
                Nb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*h^(dim-1)));
                xLim = C(1) + [-sqrt(3)*R sqrt(3)*R];
                yLim = C(2) + [-sqrt(3)*R sqrt(3)*R];
                zLim = C(3) + [-sqrt(3)*R sqrt(3)*R];
                % xB = [linspace(xLim(1),xLim(2),floor(Nb/12))', linspace(yLim(1),yLim(2),floor(Nb/12))', zlim(1);...
                %       linspace(xLim(1),xLim(2),floor(Nb/12))', linspace(yLim(1),yLim(2),floor(Nb/12))', zlim(1);...
                %       ];
            end
        end
        Nb = length(xB);
        % Interior point generation
        x = 2*(dimLCoeff(dim)*R+dimLCoeff(dim)*n^(1/dim)*h*extCoeff)*(halton(N-Nb,dim)-0.5);
        r2 = sqrt(sum(x.^2,2));
        x = x + C;
        %
        % Domain volume and area measures used for the scaling of the LS
        % problem
        %
        Vol = (dimACoeff(dim)*R^dim);
        Area = (dimPCoeff(dim)*R^(dim-1));
    end
end

% figure()
% xBplot = C + [R*cos(linspace(0,2*pi,1000))', R*sin(linspace(0,2*pi,1000))'];
% plot(xBplot(:,1),xBplot(:,2),'k-',"LineWidth",2.5)
% hold on
% plot(xc(:,1),xc(:,2),'rx','LineWidth',2,'MarkerSize',10)
% plot(xe(:,1),xe(:,2),'k.','LineWidth',2,'MarkerSize',15)
% plot(xeB(:,1),xeB(:,2),'b.','LineWidth',2,'MarkerSize',15)
% ax = gca;
% ax.FontSize = 18;
% xlabel("x","Interpreter","latex","FontSize",24)
% ylabel("y","Interpreter","latex","FontSize",24,'Rotation',0)
% legend('Domain, $$\Omega$$','Center points','Interior eval points','Boundary eval points',"Interpreter","latex",'Location','south','NumColumns',2,'Orientation','horizontal','FontSize',18)
% xlim(ax, [-1.2 1.2]);
% ylim(ax, [-1.8 1.2]);
% axis equal