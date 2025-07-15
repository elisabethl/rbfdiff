close all
clear all
setPaths;
%
% An unfitted least squares RBF-FD example for solving the Poisson equation in a dim 
% dimensional ball using the available RBF routines.
%
dim = 1;                            % dim = 1,2 or 3
display = 1;                        % Plot solution
domain = 'ball';                    % ball or cube
scaling = 1;                        % Include scaling of the LS problem
mvCentres = 1;                      % Option to have a Y point on top of all X points inside
q = 2;                              % Oversampling
h = 0.01;                            % approximate fill distance
ep = 1;                             % Not relevant for 'r3' basis
phi = 'r3';                         % Choice of basis 
pdeg = 2;                           % Polynomial extension
if pdeg == -1
    n = 2*nchoosek(5+dim,dim);
else
    n = 2*nchoosek(pdeg+dim,dim);   % Stencil size
end
extCoeff = 0;                     % Extension size (in % of stencil)
%
% Place N centre points on line/circle/sphere, start with boundary if using fitted method
%
C = zeros(1,dim);
R = 1;
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];   % Ratios between rectangle/circle, cube/sphere area and volume
dimACoeff = [1 pi (4/3)*pi];            % constant for size of domain
dimPCoeff = [1 2*pi 4*pi];              % constant for computing size of boundary
N = floor((R+n^(1/dim)*h*extCoeff)/(0.5*h))^dim;
Ncube = ceil(dimRat(dim)*(N));
xc = 2*(R+n^(1/dim)*h*extCoeff)*(halton(Ncube,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=(R+n^(1/dim)*h*extCoeff));
pos = pos(1:N);
xcBlock = xc;
xc = xc(pos,:) + C; 
%
% Place M evaluation points. Start with boundary points
%
M = ceil(N*q);
hy = ((dimACoeff(dim)*R^dim)/M).^(1/dim);
if dim == 1
    xeB = [-R, R]';
elseif dim == 2
    Mb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*hy^(dim-1)));
    xeB = [R*cos(linspace(-pi,pi,Mb)'), R*sin(linspace(-pi,pi,Mb)')];
elseif dim == 3
    Mb = ceil((dimPCoeff(dim)*R^(dim-1))/(dimACoeff(dim-1)*hy^(dim-1)));
    % Fibonacci points on sphere
    ratio = 1+sqrt(5);
    ind = [0:Mb-1]' + 0.5;
    theta = pi*ratio*ind;
    fi = acos(1-(2*ind)/(Mb));
    xeB = [R.*sin(fi).*cos(theta), R.*sin(fi).*sin(theta), R.*cos(fi)];
end
Mb = length(xeB);
Ncube = ceil(dimRat(dim)*(M-Mb));
xe = 2*R*(halton(Ncube,dim)-0.5);
r2 = sqrt(sum(xe.^2,2));
pos = find(r2<=R);
pos = pos(1:(M-Mb));
xe = xe(pos,:) + C; 
Min = length(xe);
%
% Move evaluation points on the inside to the closest centre point
%
if mvCentres
    xcIn = xc(sqrt(sum(xc.^2,2))<=R,:);
    xeTemp = xe;
    for i = 1:length(xcIn)
        idMove = knnsearch(xeTemp,xcIn(i,:),'K',1);
        xe(idMove,:) = xcIn(i,:);
        xeTemp(idMove,:) = R.*ones(1,dim);
    end
    idMove = knnsearch(xe,xcIn,'K',1);
    xe(idMove,:) = xcIn;
end
xeAll = [xe; xeB]; % Interior and boundary points (all centres for collocation)
%
% Ensure center points are not too far from boundary and make evaluation point - stencil list 
%
if extCoeff ~= 0
    xc = xc(unique(knnsearch(xc,xeAll,'K',ceil(extCoeff*n))),:);
end
N = length(xc);
ptStencilList = knnsearch(xc,xeAll,'K',1);
%
% Constructing global LS-RBF-FD approximation to evaluation, Laplace and boundary operators Min x N
%
Eglobal = spalloc(M,N,M*n);
Lglobal = spalloc(Min,N,Min*n);
LtestGlob = spalloc(M,N,M*n);
Bglobal = spalloc(Mb,N,Mb*n);
for i = 1:N
    [id,~] = knnsearch(xc,xc(i,:),'K',n);
    xcLoc = xc(id,:); % stencil
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    idInt = find(ptStencilList(1:Min)==i);
    idBnd = find(ptStencilList(Min+1:M)==i);
    % Compute evaluation, Laplace and boundary operators
    E = RBFDiffMat(0,Psi,xeAll([idInt; idBnd+Min],:));
    L = RBFDiffMat(1.5,Psi,xeAll(idInt,:));
    B = RBFDiffMat(0,Psi,xeAll(idBnd+Min,:));
    
    Eglobal([idInt; idBnd+Min],id) = E + Eglobal([idInt; idBnd+Min],id);
    Lglobal(idInt,id) = L + Lglobal(idInt,id);
    Bglobal(idBnd,id) = B + Bglobal(idBnd,id);
end
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xe(:,1)); fun(xeB(:,1))];
    uExact = fun(xeAll(:,1));
    ucExact = fun(xc(:,1));
    lapAnalytic = lapFun(xe(:,1));
    bndAnalytic = [0; 0];
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(:,1),xe(:,2)); fun(xeB(:,1),xeB(:,2))];
    uExact = fun(xeAll(:,1),xeAll(:,2));
    ucExact = fun(xc(:,1),xc(:,2));
    lapAnalytic = lapFun(xe(:,1),xe(:,2));
    bndAnalytic = fun(xeB(:,1),xeB(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xe(:,1),xe(:,2),xe(:,3)); fun(xeB(:,1),xeB(:,2),xeB(:,3))];
    uExact = fun(xeAll(:,1),xeAll(:,2),xeAll(:,3));
    ucExact = fun(xc(:,1),xc(:,2),xc(:,3));
    lapAnalytic = lapFun(xe(:,1),xe(:,2),xe(:,3));
    bndAnalytic = fun(xeB(:,1),xeB(:,2),xeB(:,3));
end 
%
% LS problem scaling
%
if scaling
    Lscale = sqrt((dimACoeff(dim)*R^dim)/Min);
    Bscale =sqrt((dimPCoeff(dim)*R^(dim-1))/Mb)*(h.^-1.5);
    A = [Lscale.*Lglobal; Bscale.*Bglobal];
    F(1:Min) = Lscale.*F(1:Min);
    F(Min+1:M) = Bscale.*F(Min+1:M);
else
    A = [Lglobal; Bglobal];
end
%
% Solution on centre points, evaluated on Y set
%
u = A\F;
ue = Eglobal*u;

if display
    plotSolution(ue,uExact,xeAll,xeB,dim);
end

%
% Operator and solution l2 errors
%
lapNumeric = Lglobal*ucExact;
evalNumeric = Bglobal*ucExact;
l2Error = norm(abs(ue-uExact),2)/norm(uExact,2);
laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
disp(['PDE error = ', num2str(l2Error)]);
disp(['Boundary Op error = ', num2str(bndError)]);
disp(['Laplace Op error = ', num2str(laplaceError)]);

%
% Plotting routines
%
function [] = plotSolution(uNumeric,uAnalytic,x,xB,dim)
    error = uNumeric-uAnalytic;
    if dim == 1
        figure()
        plot(x,uAnalytic,'b-'); 
        ax = gca;
        ax.FontSize = 18;
        xlabel("x","Interpreter","latex","FontSize",24)
        ylabel("$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
    
        figure()
        plot(x,uNumeric,'.');
        title("Numerical PDE solution")
    
        figure()
        plot(x,error,'.');
        title("Error")
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
        figure()
        scatter3(x(:,1),x(:,2),x(:,3),[],uAnalytic,'filled'); axis equal
        colorbar
        title("Exact solution")
    
        figure()
        scatter3(x(:,1),x(:,2),x(:,3),[],uNumeric,'filled'); axis equal
        colorbar
        title("Numerical solution")
    
        figure()
        scatter3(x(:,1),x(:,2),x(:,3),[],error,'filled'); axis equal
        colorbar
        title("Error")
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