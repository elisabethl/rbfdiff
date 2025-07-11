close all
clear all
setPaths;

dim = 3;
N = 5000;
q = 3;
Ne = q*N;
ep = 1;
phi = 'r3';
pdeg = 5;
n = 2*nchoosek(pdeg+dim,dim); % stencil size

C = zeros(1,dim);
R = 1;
%
% We place halton points in a circle centred at C and stretched with R
%
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
Nin = ceil(dimRat(dim)*N);
NinE = ceil(dimRat(dim)*Ne);
xc = 2*R*(halton(Nin,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=R);
pos = pos(1:N);
xc = xc(pos,:) + C; 
[~,dist] = knnsearch(xc,xc,'K',2);

if dim == 1
    xcB = [-R, R];
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

xcAll = [xc; xcB];
%
% We do the same for the evaluation points
%
% xe = 2*R*(halton(NinE,dim)-0.5);
% r2 = sqrt(sum(xe.^2,2));
% pos = find(r2<=R);
% pos = pos(1:Ne);
% xe = xe(pos,:) + C;
% [~,dist] = knnsearch(xe,xe,'K',2);
% Nb = ceil(2*pi*R/max(dist(:,2)));
% xcB = [R*cos(linspace(0,2*pi,Nb)'), R*sin(linspace(0,2*pi,Nb)')];
% xe = [xe; xcB];

Lglobal = zeros(N,length(xcAll));
for i = 1:N
    [id,~] = knnsearch(xcAll,xc(i,:),'K',n);
    xcLoc = xcAll(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    % Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc);

    L = RBFDiffMat(1.5,Psi,xcLoc(1,:));

    Lglobal(i,id) = L + Lglobal(i,id);
end
Bglobal = zeros(length(xcB),length(xcAll));
for i = 1:length(xcB)
    [id,~] = knnsearch(xcAll,xcB(i,:),'K',n);
    xcLoc = xcAll(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    % Psi = RBFInterpMat(phi,pdeg,ep,xcLox,xcLoc);

    B = RBFDiffMat(0,Psi,xcLoc(1,:));

    Bglobal(i,id) = B + Bglobal(i,id);
end

Eglobal = zeros(length(xcAll),length(xcAll));
for i = 1:length(xcAll)
    [id,~] = knnsearch(xcAll,xcAll(i,:),'K',n);
    xcLoc = xcAll(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    % Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc);

    E = RBFDiffMat(0,Psi,xcLoc(1,:));

    Eglobal(i,id) = E + Eglobal(i,id);

end

Ltest = zeros(length(xc),length(xc));
for i = 1:length(xc)
    [id,~] = knnsearch(xc,xc(i,:),'K',n);
    xcLoc = xc(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc(1,:),max(sqrt(sum((xcLoc-xcLoc(1,:)).^2,2))));
    % Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc);

    L = RBFDiffMat(1.5,Psi,xcLoc(1,:));

    Ltest(i,id) = L + Ltest(i,id);

end

if dim == 1
    fun = @(x) sin(2.*pi.*x);
    lapFun = @(x) -4.*(x.^2).*(pi.^2).*sin(2.*pi.*x);
    F = [lapFun(xc(:,1)); fun(xcB(:,1))];
    uc = fun(xcAll(:,1));
    lapAnalytic = lapFun(xc(:,1));
elseif dim == 2
    fun = @(x,y) sin(2.*pi.*x.*y);
    lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);
    F = [lapFun(xc(:,1),xc(:,2)); fun(xcB(:,1),xcB(:,2))];
    uc = fun(xcAll(:,1),xcAll(:,2));
    lapAnalytic = lapFun(xc(:,1),xc(:,2));
elseif dim == 3
    fun = @(x,y,z) sin(2.*pi.*x.*y.*z);
    lapFun = @(x,y,z) - 4.*(x.^2).*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z) - 4.*(y.^2).*(z.^2).*(pi.^2).*sin(2.*pi.*x.*y.*z);
    F = [lapFun(xc(:,1),xc(:,2),xc(:,3)); fun(xcB(:,1),xcB(:,2),xcB(:,3))];
    uc = fun(xcAll(:,1),xcAll(:,2),xcAll(:,3));
    lapAnalytic = lapFun(xc(:,1),xc(:,2),xc(:,3));
end 

A = [Lglobal; Bglobal];

u = A\F;
lapNumeric = Ltest*uc(1:N);
uNumeric = Eglobal*uc;

if dim == 1
elseif dim == 2
    figure()
    scatter(xcAll(:,1),xcAll(:,2),[],uc,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter(xcAll(:,1),xcAll(:,2),[],u,'filled'); axis equal
    colorbar
    title("Numeric PDE solution")
    
    figure
    scatter(xc(:,1),xc(:,2),[],lapAnalytic,'filled'); axis equal
    colorbar
    title("Analytical laplace")

    figure()
    scatter(xc(:,1),xc(:,2),[],lapNumeric,'filled'); axis equal
    colorbar
    title("Numeric laplace")
    
    figure()
    scatter(xcAll(:,1),xcAll(:,2),[],uNumeric,'filled'); axis equal
    colorbar
    title("Interpolated solution")
elseif dim == 3
    figure()
    scatter3(xcAll(:,1),xcAll(:,2),xcAll(:,3),[],uc,'filled'); axis equal
    colorbar
    title("Exact solution")

    figure()
    scatter3(xcAll(:,1),xcAll(:,2),xcAll(:,3),[],u,'filled'); axis equal
    colorbar
    title("Numeric PDE solution")
    
    figure
    scatter3(xc(:,1),xc(:,2),xc(:,3),[],lapAnalytic,'filled'); axis equal
    colorbar
    title("Analytical laplace")
    
    figure()
    scatter3(xc(:,1),xc(:,2),xc(:,3),[],lapNumeric,'filled'); axis equal
    colorbar
    title("Numeric laplace")
    
    figure()
    scatter3(xcAll(:,1),xcAll(:,2),xcAll(:,3),[],uNumeric,'filled'); axis equal
    colorbar
    title("Interpolated solution")
end

error = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
disp(['Laplacian error = ', num2str(error)])

error = norm(uc-u,2)/norm(uc,2);
disp(['PDE error = ', num2str(error)]);

error = norm(uc-uNumeric,2)/norm(uc,2);
disp(['Interpolation error = ', num2str(error)]);