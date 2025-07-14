addpath oldqrcode
% Compare the differentiation accuracy between the old and new implementations
clear
close all
%
% Parameters to choose
%
ep = 0.1;
dim = 3;
N = 56;
Ne = 200;
ttype = 1;
if (ttype==1)
    r = logspace(-16,-5,Ne)'; % For errors near r=0
    rplot = r;
    rstr = 'r';
elseif (ttype==2)
    r = 1-logspace(-16,-5,Ne)'; % For errors near r=1
    rplot = 1-r;
    rstr = 'Distance from $r=1$';
end

fnum = 1;

if (fnum==1)
    fun = @(x) ones(size(x,1),1);
    funx = @(x) zeros(size(x,1),1);
    funxx = @(x) zeros(size(x,1),1);
end    
%
% We let the domain be a disc or a sphere centered in the origin 
% The node points are Halton nodes without clustering
%
C = zeros(1,dim);
R = 1;
dimRat = 1.2*8/(4*pi/3);
Nin = ceil(dimRat*N);
xc = 2*(halton(Nin,dim)-0.5);
r2 = sum(xc.^2,2);
pos = find(r2<=R);
pos = pos(1:N);
xc = xc(pos,:); 
%
% Evaluate on a ray for some phi and theta, such that we can plot the error
%
phi = 0.234*pi*ones(Ne,1);
theta = 0.9*pi*ones(Ne,1);

[xe(:,1),xe(:,2),xe(:,3)] = sph2cart(theta,phi,r);

uc = fun(xc);
ue = fun(xe);
uxe = funx(xe);
uxxe = funxx(xe);
%
% Initialize the approximation
%
Psi = RBF_QR_diffmat_3D(0,xc,ep,C,R);
[~,PsiO] = RBF_QR_diffmat_3D_old('0',xe,xc,ep);
%
% Compute the differentiation matrices
%
[B,Te] = RBF_QR_diffmat_3D(0,Psi,xe);
[Bx,Te] = RBF_QR_diffmat_3D(1,Psi,Te);
[BL,Te] = RBF_QR_diffmat_3D(1.5,Psi,Te);
%
BO = RBF_QR_diffmat_3D_old('0',xe,PsiO);
BxO = RBF_QR_diffmat_3D_old('x',xe,PsiO);
ByO = RBF_QR_diffmat_3D_old('y',xe,PsiO);
BLO = RBF_QR_diffmat_3D_old('L',xe,PsiO);
%
% Apply them to the right hand side
%
u = B*uc;
uO = BO*uc;
Lu = BL*uc;
LuO = BLO*uc;
%
% Compute errors
%
erru = u-ue;
erruO = uO-ue;
errL = Lu-uxxe; % In this case
errLO = LuO-uxxe; % In this case


figure
H = loglog(rplot,abs(erru),'-',rplot,abs(erruO),'--');
set(H,'LineWidth',1.5)
adjustFigProp
xlabel(rstr,'interpreter','latex')
title('Interpolation error','Interpreter','latex')

figure
H = loglog(rplot,abs(errL),'-',rplot,abs(errLO),'--');
set(H,'LineWidth',1.5)
adjustFigProp
xlabel(rstr,'interpreter','latex')
title(['Error in Laplacian, $N=$' num2str(N) ', $\varepsilon=$' num2str(ep)],'Interpreter','latex')

% In title, we want ep, 
function adjustFigProp()
set(gca,'LineWidth',1.5,'FontSize',18)
set(gcf,'Color','w')
axis tight
grid
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
set(gca,'XTick',[1e-15 1e-10 1e-5])
set(gca,'YTick',[1e-8 1e-6 1e-4 1e-2 1e0 1e2 1e4])
legend('New code','Old code')
end


