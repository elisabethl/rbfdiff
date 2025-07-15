close all
clear all
setPaths;
%
% A collocation RBF-FD example for solving the Poisson equation in a dim 
% dimensional ball using the available RBF routines.
%
dim = 2;                        % dim = 1,2 or 3
h = 0.1;                      % approximate fill distance
ep = 1;                         % Not relevant for 'r3' basis
phi = 'r3';
pdeg = 4;
n = 2*nchoosek(pdeg+dim,dim);   % stencil size
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
% Constructing global RBF-FD approximation to Laplace operator Nin x N
%
Lglobal = spalloc(Nin,N,Nin*n);
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
Bglobal = spalloc(Nb,N,Nb*n);
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
    ax = gca;
    ax.FontSize = 18;
    xlabel("x","Interpreter","latex","FontSize",24)
    ylabel("y","Interpreter","latex","FontSize",24)
    zlabel("$$u_E$$","Interpreter","latex","FontSize",24,'Rotation',0)
    set(G,'EdgeColor','none')
    shading interp

    figure()
    G=trisurf(T,xAll(:,1),xAll(:,2),u);
    hold on
    plot(xcB(:,1),xcB(:,2),'k-',"LineWidth",1.5)
    ax = gca;
    ax.FontSize = 18;
    xlabel("x","Interpreter","latex","FontSize",24)
    ylabel("y","Interpreter","latex","FontSize",24)
    zlabel("$$u$$","Interpreter","latex","FontSize",24,'Rotation',0)
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

    figure()
    xBplot = C + [R*cos(linspace(0,2*pi,1000))', R*sin(linspace(0,2*pi,1000))'];
    plot(xBplot(:,1),xBplot(:,2),'k-',"LineWidth",2.5)
    hold on
    plot(xc(:,1),xc(:,2),'rx','LineWidth',2,'MarkerSize',10)
    plot(xcB(:,1),xcB(:,2),'b.','LineWidth',2,'MarkerSize',15)
    ax = gca;
    ax.FontSize = 18;
    xlabel("x","Interpreter","latex","FontSize",24)
    ylabel("y","Interpreter","latex","FontSize",24,'Rotation',0)
    legend('Domain, $$\Omega$$','Center points','Boundary points',"Interpreter","latex",'Location','south','Orientation','horizontal','FontSize',18)
    xlim(ax, [-1.4 1.2]);
    ylim(ax, [-1.4 1.2]);
    axis equal
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