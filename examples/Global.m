close all
clear all
setPaths;
%
% A collocation RBF-FD example for solving the Poisson equation in a dim 
% dimensional ball using the available RBF routines.
%
dim = 2;                        % dim = 1,2 or 3
h = 0.1;                        % approximate fill distance
ep = 0.1;                         % Not relevant for 'r3' basis
phi = 'rbfqr';
pdeg = -1;
%
% We place halton points in a line/circle/sphere centred at C with radius R
%
C = zeros(1,dim);
R = 1;
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
Nin = floor(1/(0.5*h))^dim;
Ncube = ceil(dimRat(dim)*Nin);
xc = 2*R*(halton(Ncube,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=R);
pos = pos(1:Nin);
xc = xc(pos,:) + C; 
[~,dist] = knnsearch(xc,xc,'K',2);
%
% Place boundary points using fill distance measure from interior
%
if dim == 1
    xcB = [-R, R]';
elseif dim == 2
    Nb = ceil(2*pi*R/max(dist(:,2)));
    xcB = [R*cos(linspace(0,2*pi,Nb)'), R*sin(linspace(0,2*pi,Nb)')];
elseif dim == 3
    % Fibonacci points on sphere
    Nb = ceil(4*pi*R/max(dist(:,2)));
    ratio = 1+sqrt(5);
    ind = [0:Nb-1]' + 0.5;
    theta = pi*ratio*ind;
    fi = acos(1-(2*ind)/(Nb));
    xcB = [R.*sin(fi).*cos(theta), R.*sin(fi).*sin(theta), R.*cos(fi)];
end

xAll = [xc; xcB]; % Interior and boundary points (all centres for collocation)
N = length(xAll);
Nb = length(xcB);
%
% Constructing global Kansa approximation to Laplace and boundary operators
%
Psi = RBFInterpMat(phi,pdeg,ep,xAll,C,R);
Lglobal = RBFDiffMat(1.5,Psi,xc);
Bglobal = RBFDiffMat(0,Psi,xcB);
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xc(:,1)); fun(xcB(:,1))];
    uc = fun(xAll(:,1));
    lapAnalytic = lapFun(xc(:,1)); 
    bndAnalytic = fun(xcB(:,1));
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xc(:,1),xc(:,2)); fun(xcB(:,1),xcB(:,2))];
    uc = fun(xAll(:,1),xAll(:,2));
    lapAnalytic = lapFun(xc(:,1),xc(:,2));  
    bndAnalytic = fun(xcB(:,1),xcB(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xc(:,1),xc(:,2),xc(:,3)); fun(xcB(:,1),xcB(:,2),xcB(:,3))];
    uc = fun(xAll(:,1),xAll(:,2),xAll(:,3));
    lapAnalytic = lapFun(xc(:,1),xc(:,2),xc(:,3));    
    bndAnalytic = fun(xcB(:,1),xcB(:,2),xcB(:,3));
end 

A = [Lglobal; Bglobal];

u = A\F;

error = u-uc;

if dim == 1
    figure()
    plot(xAll,uc,'.'); 
    title("Exact solution")

    figure()
    plot(xAll,u,'.');
    title("Numerical PDE solution")

    figure()
    plot(xAll,error,'.');
    title("Error")
elseif dim == 2
    T = delaunay(xAll(:,1),xAll(:,2));
    G=trisurf(T,xAll(:,1),xAll(:,2),uc);
    hold on
    plot(xcB(:,1),xcB(:,2),'k-',"LineWidth",1.5)
    G=trisurf(T,xAll(:,1),xAll(:,2),uc);
    set(G,'EdgeColor','none')
    shading interp

    figure()
    G=trisurf(T,xAll(:,1),xAll(:,2),u);
    hold on
    plot(xcB(:,1),xcB(:,2),'k-',"LineWidth",1.5)
    set(G,'EdgeColor','none')
    shading interp

    figure()
    G=trisurf(T,xAll(:,1),xAll(:,2),error);
    hold on
    plot(xcB(:,1),xcB(:,2),'k-',"LineWidth",1.5)
    ax = gca;
    ax.FontSize = 18;
    xlabel("x","Interpreter","latex","FontSize",24)
    ylabel("y","Interpreter","latex","FontSize",24)
    zlabel("$$u-u_E$$","Interpreter","latex","FontSize",24)
    set(G,'EdgeColor','none')
    shading interp
elseif dim == 3
    figure()
    scatter3(xAll(:,1),xAll(:,2),xAll(:,3),[],uc,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter3(xAll(:,1),xAll(:,2),xAll(:,3),[],u,'filled'); axis equal
    colorbar
    title("Numerical solution")

    figure()
    scatter3(xAll(:,1),xAll(:,2),xAll(:,3),[],error,'filled'); axis equal
    colorbar
    title("Error")
end
%
% Operator and solution l2 errors
%
lapNumeric = Lglobal*uc;
evalNumeric = Bglobal*uc;
l2Error = norm(abs(error),2)/norm(uc,2);
laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
bndError = norm(evalNumeric-bndAnalytic,2)/(norm(bndAnalytic,2) + double(max(abs(bndAnalytic))==0));
disp(['PDE error = ', num2str(l2Error)]);
disp(['Boundary Op error = ', num2str(bndError)]);
disp(['Laplace Op error = ', num2str(laplaceError)]);