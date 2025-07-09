close all
clear all
setPaths;

dim = 2;
N = 100;
q = 3;
Ne = q*N;
n = 5; % stencil size
ep = 0.1;
phi = 'r3';
pdeg = 2;

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

for i = 1:N
    [id,~] = knnsearch(xc,xc(1,:),'K',n);
    xcLoc = xc(id,:);
    Psi = RBFInterpMat(phi,pdeg,ep,xcLoc);
end
