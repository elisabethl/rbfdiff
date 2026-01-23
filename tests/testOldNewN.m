addpath tests/oldqrcode
%
% This is a test to see if there is a difference in the result from the old
% and new code, but with a specific focus on changing N to extreme sizes.
%
clear
Nvec = logspace(log10(15),log10(4000),50); % The number of node points in the domain
ep = 1e-1; % A fixed small shape parameter

fnum = 6;

% Add some functions to choose from and compute their derivatives symbolically
syms x y
if fnum==1
    f = 0.5*(x^0+y^0);
elseif fnum==2
    f = x + y;
elseif fnum==3
    f = x^2 + y^2;
elseif fnum==4
    f = sin(x-2*y);
elseif fnum==5
    f = 25/(1+25*(x^2+0.5*y^2));
elseif fnum==6
    f = sin(x*y+y^2)-cos(x^2);
end   
fun = matlabFunction(f,'vars',{'x','y'});
funx = matlabFunction(diff(f,x),'vars',{'x','y'});
funy = matlabFunction(diff(f,y),'vars',{'x','y'});
funxx = matlabFunction(diff(f,x,2),'vars',{'x','y'});
funxy = matlabFunction(diff(f,x,y),'vars',{'x','y'});
funyy = matlabFunction(diff(f,y,2),'vars',{'x','y'});

for j=1:length(Nvec)
    N = round(Nvec(j))
    Ne = min(400,2*N); % The number of evaluation points
    clear xk xe
    % Use Halton nodes in the domain, clustered towards the boundary
    d=2;
    Nin = ceil(1.2*4/pi*N);
    xk = 2*(halton(Nin,d)-0.5);
    [th,r]=cart2pol(xk(:,1),xk(:,2));
    pos = find(r<=1);
    pos = pos(1:N);
    xk = xk(pos,:); th = th(pos); r = r(pos);
    r = sin(pi/2*r);
    [xk(:,1),xk(:,2)] = pol2cart(th,r);

    % Evaluation points (Ne = nt*nr = 3*n*n)
    nr = ceil(sqrt(Ne/3));
    nt = 3*nr;
    th = 0:2*pi/nt:2*pi; % First and last is the same
    r = 0:1/(nr-1):1;
    [tt,rr] = meshgrid(th,r);
    sz = size(tt);
    [xe(:,1), xe(:,2)] = pol2cart(tt(:),rr(:));
    %
    % Apply to a function and compare with exact solution
    %
    uk = fun(xk(:,1),xk(:,2));
    if length(uk)==1
        uk = uk*ones(size(xk(:,1)));
    end    
    
    [A,Psi,Te] = RBF_QR_diffmat_2D(0,xk,ep,xe);
    [G,Te] = RBF_QR_diffmat_2D(1,Psi,Te);
    [L,Te] = RBF_QR_diffmat_2D(1.5,Psi,Te);
    [H,Te] = RBF_QR_diffmat_2D(2,Psi,Te);
    %
    % Do the same with the old code
    %
    [A2,Psi]=RBF_QR_diffmat_2D_old('0',xe,xk,ep);
    [G2{1},Psi]=RBF_QR_diffmat_2D_old('x',xe,Psi);
    [G2{2},Psi]=RBF_QR_diffmat_2D_old('y',xe,Psi);
    [L2,Psi]=RBF_QR_diffmat_2D_old('L',xe,Psi);
    [H2{1,1},Psi]=RBF_QR_diffmat_2D_old('xx',xe,Psi);
    [H2{1,2},Psi]=RBF_QR_diffmat_2D_old('xy',xe,Psi);
    [H2{2,2},Psi]=RBF_QR_diffmat_2D_old('yy',xe,Psi);
    fe = fun(xe(:,1),xe(:,2));
    fxe = funx(xe(:,1),xe(:,2));
    fye = funy(xe(:,1),xe(:,2));
    fxxe = funxx(xe(:,1),xe(:,2));
    fxye = funxy(xe(:,1),xe(:,2));
    fyye = funyy(xe(:,1),xe(:,2));
    fLe = funxx(xe(:,1),xe(:,2))+funyy(xe(:,1),xe(:,2));
    maxerr(j) = max(abs(A*uk-fe))./max(abs(fe));
    errx(j) = max(abs(G{1}*uk-fxe))./max(abs(fxe));
    erry(j) = max(abs(G{2}*uk-fye))./max(abs(fye));
    errxx(j) = max(abs(H{1,1}*uk-fxxe))./max(abs(fxxe));
    errxy(j) = max(abs(H{1,2}*uk-fxye))./max(abs(fxye));
    erryy(j) = max(abs(H{2,2}*uk-fyye))./max(abs(fyye));
    errL(j) = max(abs(L*uk-fLe))./max(abs(fLe));
    %
    % Old errors
    %
    maxerr2(j) = max(abs(A2*uk-fe))./max(abs(fe));
    errx2(j) = max(abs(G2{1}*uk-fxe))./max(abs(fxe));
    erry2(j) = max(abs(G2{2}*uk-fye))./max(abs(fye));
    errxx2(j) = max(abs(H2{1,1}*uk-fxxe))./max(abs(fxxe));
    errxy2(j) = max(abs(H2{1,2}*uk-fxye))./max(abs(fxye));
    erryy2(j) = max(abs(H2{2,2}*uk-fyye))./max(abs(fyye));
    errL2(j) = max(abs(L2*uk-fLe))./max(abs(fLe));
end

figure
H1=loglog(Nvec,maxerr,Nvec,errx,Nvec,erry,Nvec,errxx,Nvec,errxy,Nvec,erryy,Nvec,errL);
hold on
H2=loglog(Nvec,maxerr2,Nvec,errx2,Nvec,erry2,Nvec,errxx2,Nvec,errxy2,Nvec,erryy2,Nvec,errL2);

set(gca,'FontSize',18)
set(gcf,'Color','w')

for k=1:length(H1)
    set(H1(k),'LineStyle','-','LineWidth',1.5)
    col=get(H1(k),'Color')
    set(H2(k),'LineStyle','none','LineWidth',1.5,'Marker','.', ...
       'MarkerSize',10,'Color',col)
end    
H=legend([H1(1) H2(1)],'New code','Old code');
set(H,'Location','NorthWest')
xlabel('The number of node points $n$','interpreter','latex')
ylabel('Maximum error','interpreter','latex')
title('Shape parameter $\varepsilon=0.1$','interpreter','latex')
grid
set(gca,'XMinorGrid','off')
set(gca,'YTick',[1e-10 1e-5 1])
