close all
clear all
%
% Testing inHouseLegendre against legendre_modified
%
addpath("rbfqr3D")
%
% x close to 0, 1, -1 or just in range [-1,1]
%
x = [-1, logspace(-16,-10,100) - 1];
% x = [1, -logspace(-16,-10,100) + 1];
% x = [0, logspace(-16,-10,50)];
% x = [-logspace(-16,-10,50), 0];
% x = linspace(-1,1,100);
%
% Associated Lagrange polynomial degree
%
l = 3;
[P, dPdth, pOverS, dPdth2, A, B] = normalizedLegendre(l, x);
[P_mod, dPdth_mod, pOverS_mod, dPdth2_mod, A_mod, B_mod] = legendre_modified(l, x, "norm");

% Error in P
for i = 1:l+1
    figure()
    loglog(x,P(i,:),'b-',...
           x,P_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function $$P^{',num2str(i-1),'}_{',num2str(l),'}$$'],'Interpreter','latex');
    error = norm(P(i,:)-P_mod(i,:),2)/norm(P_mod(i,:),2);
    disp(['l2 error in P, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end

% Error in dPdth
for i = 1:l+1
    figure()
    loglog(x,dPdth(i,:),'b-',...
           x,dPdth_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function derivative  $$dP^{',num2str(i-1),'}_{',num2str(l),'}/d\theta$$'],'Interpreter','latex');
    error = norm(dPdth(i,:)-dPdth_mod(i,:),2)/norm(dPdth_mod(i,:),2);
    disp(['l2 error in dPdth, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end

% Error in pOverS
for i = 1:l+1
    figure()
    loglog(x,pOverS(i,:),'b-',...
           x,pOverS_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function $$P^{',num2str(i-1),'}_{',num2str(l),'}/\sqrt{1-x^2}$$'],'Interpreter','latex');
    error = norm(pOverS(i,:)-pOverS_mod(i,:),2)/norm(pOverS_mod(i,:),2);
    disp(['l2 error in pOverS, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end

% Error in dPdth2
for i = 1:l+1
    figure()
    loglog(x,dPdth2(i,:),'b-',...
           x,dPdth2_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function derivative $$d^2P^{',num2str(i-1),'}_{',num2str(l),'}/d\theta^2$$'],'Interpreter','latex');
    error = norm(dPdth2(i,:)-dPdth2_mod(i,:),2)/norm(dPdth2_mod(i,:),2);
    disp(['l2 error in dPdth2, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end

% Error in A
for i = 1:l+1
    figure()
    loglog(x,A(i,:),'b-',...
           x,A_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function combination $$A$$'],'Interpreter','latex');
    error = norm(A(i,:)-A_mod(i,:),2)/norm(dPdth2_mod(i,:),2);
    disp(['l2 error in A, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end

% Error in B
for i = 1:l+1
    figure()
    loglog(x,B(i,:),'b-',...
           x,B_mod(i,:),'r--')
    xlim([min(x),max(x)])
    title(['Associated Legendre function combination $$B$$'],'Interpreter','latex');
    error = norm(dPdth2(i,:)-dPdth2_mod(i,:),2)/norm(dPdth2_mod(i,:),2);
    disp(['l2 error in B, for ','m = ',num2str(i-1),' is ',num2str(error)]);
end