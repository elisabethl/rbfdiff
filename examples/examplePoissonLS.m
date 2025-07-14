close all
clear all
setPaths;
%
% A least squares, unfitted RBF-FD example for solving the Poisson equation in a dim 
% dimensional ball using the available RBF routines.
%
dim = 2;                        % dim = 1,2 or 3
scaling = 1;                    % Include scaling of the LS problem
q = 2;                          % Oversampling
h = 0.1;                        % approximate fill distance
ep = 1;                         % Not relevant for 'r3' basis
phi = 'r3';
pdeg = 2;
n = 2*nchoosek(pdeg+dim,dim);   % stencil size
extCoeff = 0.5;
%
% We place centre points in a line/circle/sphere extended by half a stencil
%
C = zeros(1,dim);
R = 1;
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
dimACoeff = [1 pi (4/3)*pi];
dimPCoeff = [1 2*pi 4*pi];
N = floor((R+n^(1/dim)*h*extCoeff)/(0.5*h))^dim;
M = N*q;
Ncube = ceil(dimRat(dim)*N);
xc = 2*(R+n^(1/dim)*h*extCoeff)*(halton(Ncube,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=(R+n*h*extCoeff));
pos = pos(1:N);
xcBlock = xc;
xc = xc(pos,:) + C; 
%
% Place evaluation points. Start with boundary points
%
if dim == 1
    xeB = [-R, R]';
elseif dim == 2
    approx_hy = ((pi*R^2)/(N*q)).^(1/dim);
    Nb = ceil(2*pi*R/approx_hy);
    xeB = [R*cos(linspace(0,2*pi,Nb)'), R*sin(linspace(0,2*pi,Nb)')];
elseif dim == 3
    % Fibonacci points on sphere
    approx_hy = ((pi*R^2)/(N*q)).^(1/dim);
    Nb = ceil(4*pi*R/approx_hy);
    ratio = 1+sqrt(5);
    ind = [0:Nb-1]' + 0.5;
    theta = pi*ratio*ind;
    fi = acos(1-(2*ind)/(Nb));
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
xcIn = xc(sqrt(sum(xc.^2,2))<=R,:);
xeTemp = xe;
for i = 1:length(xcIn)
    idMove = knnsearch(xeTemp,xcIn(i,:),'K',1);
    xe(idMove,:) = xcIn(i,:);
    xeTemp(idMove,:) = R.*ones(1,dim);
end
idMove = knnsearch(xe,xcIn,'K',1);
xe(idMove,:) = xcIn;
xeAll = [xe; xeB]; % Interior and boundary points (all centres for collocation)
%
% Ensure center points are not too far from boundary and make evaluation point - stencil list 
%
xc = xc(unique(knnsearch(xc,xeAll,'K',ceil(extCoeff*n))),:);
ptStencilList = knnsearch(xc,xeAll,'K',1);
N = length(xc);
%
% Constructing global LS-RBF-FD approximation to evaluation, Laplace and boundary operators Min x N
%
Eglobal = zeros(M,N);
Lglobal = zeros(Min,N);
LtestGlob = zeros(M,N);
Bglobal = zeros(Mb,N);
for i = 1:N
    [id,~] = knnsearch(xc,xc(i,:),'K',n);
    xcLoc = xc(id,:); % stencil
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    idInt = find(ptStencilList(1:Min)==i);
    idBnd = find(ptStencilList(Min+1:M)==i);
    E = RBFDiffMat(0,Psi,xeAll([idInt; idBnd+Min],:));
    Ltest = RBFDiffMat(1.5,Psi,xeAll([idInt; idBnd+Min],:));
    L = RBFDiffMat(1.5,Psi,xeAll(idInt,:));
    B = RBFDiffMat(0,Psi,xeAll(idBnd+Min,:));
    
    Eglobal([idInt; idBnd+Min],id) = E + Eglobal([idInt; idBnd+Min],id);
    Lglobal(idInt,id) = L + Lglobal(idInt,id);
    LtestGlob([idInt; idBnd+Min],id) = Ltest + LtestGlob([idInt; idBnd+Min],id);
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
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xe(:,1),xe(:,2)); fun(xeB(:,1),xeB(:,2))];
    uExact = fun(xeAll(:,1),xeAll(:,2));
    ucExact = fun(xc(:,1),xc(:,2));
    lapAnalytic = lapFun(xe(:,1),xe(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(x.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xe(:,1),xe(:,2),xe(:,3)); fun(xeB(:,1),xeB(:,2),xeB(:,3))];
    uExact = fun(xeAll(:,1),xeAll(:,2),xeAll(:,3));
    ucExact = fun(xc(:,1),xc(:,2),xc(:,3));
    lapAnalytic = lapFun(xe(:,1),xe(:,2),xe(:,3));
end 

if scaling == 1
    Lscale = sqrt((dimACoeff(dim)*R^dim)/Min);
    Bscale =sqrt((dimPCoeff(dim)*R^(dim-1))/Mb)*(h.^-1.5);
    % Lscale = sqrt(1/Min);
    % Bscale = sqrt(1/Mb)*(h.^-1);
    A = [Lscale.*Lglobal; Bscale.*Bglobal];
    F(1:Min) = Lscale.*F(1:Min);
    F(Min+1:M) = Bscale.*F(Min+1:M);
elseif scaling == 2
    Lscale = sqrt(sum(Lglobal.^2,2));
    Bscale = sqrt(sum(Bglobal.^2,2));
    A = [Lscale.*Lglobal; Bscale.*Bglobal];
    F(1:Min) = Lscale.*F(1:Min);
    F(Min+1:M) = Bscale.*F(Min+1:M);
else
    A = [Lglobal; Bglobal];
end

u = A\F;
ue = Eglobal*u;

error = abs(ue-uExact);

if dim == 1
    figure()
    plot(xeAll,uExact,'.'); 
    title("Exact solution")

    figure()
    plot(xeAll,ue,'.');
    title("Numerical PDE solution")

    figure()
    plot(xeAll,error,'.');
    title("Error")
elseif dim == 2
    figure()
    scatter(xeAll(:,1),xeAll(:,2),[],uExact,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter(xeAll(:,1),xeAll(:,2),[],ue,'filled'); axis equal
    colorbar
    title("Numerical PDE solution")

    figure()
    scatter(xeAll(:,1),xeAll(:,2),[],error,'filled'); axis equal
    colorbar
    hold on
    plot(xc(:,1),xc(:,2),'rx')
    title("Error")
elseif dim == 3
    figure()
    scatter3(xeAll(:,1),xeAll(:,2),xeAll(:,3),[],uExact,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter3(xeAll(:,1),xeAll(:,2),xeAll(:,3),[],ue,'filled'); axis equal
    colorbar
    title("Numerical solution")

    figure()
    scatter3(xeAll(:,1),xeAll(:,2),xeAll(:,3),[],error,'filled'); axis equal
    colorbar
    title("Error")
end

lapNumeric = Lglobal*ucExact;
evalNumeric = Bglobal*ucExact;
l2Error = norm(error,2)/norm(uExact,2);
laplaceError = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
bndError = norm(evalNumeric-fun(xeB(:,1),xeB(:,2)),2)/norm(fun(xeB(:,1),xeB(:,2)),2);
disp(['PDE error = ', num2str(l2Error)]);
disp(['Boundary error = ', num2str(bndError)]);
disp(['Laplace error = ', num2str(laplaceError)]);