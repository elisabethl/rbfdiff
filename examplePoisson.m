close all
clear all
setPaths;

dim = 2;
N = 500;
q = 3;
Ne = q*N;
ep = 1;
phi = 'r3';
pdeg = 4;
n = 2*nchoosek(pdeg+dim,dim); % stencil size

C = zeros(1,dim);
R = 1;
%
% We place halton points in a circle centred at "shift" and stretched with R
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
Nb = ceil(2*pi*R/max(dist(:,2)));
xcB = [R*cos(linspace(0,2*pi,Nb)'), R*sin(linspace(0,2*pi,Nb)')];
xc = [xc; xcB];
N = length(xc);
%
% We do the same for the evaluation points
%
xe = 2*R*(halton(NinE,dim)-0.5);
r2 = sqrt(sum(xe.^2,2));
pos = find(r2<=R);
pos = pos(1:Ne);
xe = xe(pos,:) + C;
[~,dist] = knnsearch(xe,xe,'K',2);
Nb = ceil(2*pi*R/max(dist(:,2)));
xcB = [R*cos(linspace(0,2*pi,Nb)'), R*sin(linspace(0,2*pi,Nb)')];
xe = [xe; xcB];

Lglobal = zeros(N,N);
for i = 1:N
    [id,~] = knnsearch(xc,xc(i,:),'K',n);
    xcLoc = xc(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc,xcLoc);

    L = RBFDiffMat(1.5,Psi,xcLoc(1));
    
    Lglobal(id(1),id) = L + Lglobal(id(1),id);
end

fun = @(x,y) sin(2.*pi.*x.*y);
lapFun = @(x,y) - 4.*(x.^2).*(pi.^2).*sin(2.*pi.*x.*y) - 4.*(y.^2).*(pi.^2).*sin(2.*pi.*x.*y);

uc = fun(xc(:,1),xc(:,2));
lapAnalytic = lapFun(xc(:,1),xc(:,2));

lapNumeric = Lglobal*uc;

figure
scatter(xc(:,1),xc(:,2),[],lapAnalytic,'filled'); axis equal
colorbar
title("Exact laplacian")
figure()
scatter(xc(:,1),xc(:,2),[],lapNumeric,'filled'); axis equal
colorbar
title("Numeric laplacian")

error = norm(lapAnalytic-lapNumeric,2)/norm(lapAnalytic,2);
disp(['Laplacian error = ', num2str(error)])