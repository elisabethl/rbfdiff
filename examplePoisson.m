close all
clear all
setPaths;
%
% A collocation RBF-FD example for solving the Poisson equation in a dim 
% dimensional ball using the available RBF routines.
%
dim = 2;                        % dim = 1,2 or 3
Nin = 1000; 
ep = 1;                         % Not relevant for 'r3' basis
phi = 'r3';
pdeg = 3;
n = 2*nchoosek(pdeg+dim,dim);   % stencil size
%
% We place halton points in a line/circle/sphere centred at C and stretched with R
%
C = zeros(1,dim);
R = 1;
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
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
% Constructing global RBF-FD approximation to Laplace operator Nin x N
%
Lglobal = zeros(Nin,N);
for i = 1:Nin
    [id,~] = knnsearch(xAll,xc(i,:),'K',n);
    xcLoc = xAll(id,:); % stencil
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));

    L = RBFDiffMat(1.5,Psi,xcLoc(1,:));

    Lglobal(i,id) = L + Lglobal(i,id);
end
%
% Constructing global RBF-FD approximation to Dirichlet boundary operator Nb x N
%
Bglobal = zeros(Nb,N);
for i = 1:Nb
    [id,~] = knnsearch(xAll,xcB(i,:),'K',n);
    xcLoc = xAll(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));

    B = RBFDiffMat(0,Psi,xcLoc(1,:));

    Bglobal(i,id) = B + Bglobal(i,id);
end
%
% Manufactured solution to construct forcing and BC
%
if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xc(:,1)); fun(xcB(:,1))];
    uc = fun(xAll(:,1));
    lapAnalytic = lapFun(xc(:,1));
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xc(:,1),xc(:,2)); fun(xcB(:,1),xcB(:,2))];
    uc = fun(xAll(:,1),xAll(:,2));
    lapAnalytic = lapFun(xc(:,1),xc(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xc(:,1),xc(:,2),xc(:,3)); fun(xcB(:,1),xcB(:,2),xcB(:,3))];
    uc = fun(xAll(:,1),xAll(:,2),xAll(:,3));
    lapAnalytic = lapFun(xc(:,1),xc(:,2),xc(:,3));
end 

A = [Lglobal; Bglobal];

u = A\F;

error = abs(u-uc);

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
    figure()
    scatter(xAll(:,1),xAll(:,2),[],uc,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter(xAll(:,1),xAll(:,2),[],u,'filled'); axis equal
    colorbar
    title("Numerical PDE solution")

    figure()
    scatter(xAll(:,1),xAll(:,2),[],error,'filled'); axis equal
    colorbar
    title("Error")
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
l2Error = norm(error,2)/norm(uc,2);
disp(['PDE error = ', num2str(l2Error)]);