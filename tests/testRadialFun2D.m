close all
%
% We test some of the radial functions needed for 2D and 3D in RBF-QR here
%
ttype = '1'; % Alternatives are 0 and one for points near these regions.
n = 500;
if (ttype=='0')
    x = logspace(-16,0,n)'; % Positive half line
else
    x = logspace(-16,0,n)'; % Clustered around 1 instead
    x = 1-x;
end    
x2 = x.^2;

jmax = 50;
jall = 0:jmax;
%
% The initial (old) implementation
%
tol = 10*eps; 
T_old = cos(acos(x(:,1))*(0:jmax));
fac  = 1-x2;     
pos = find(abs(x-1)<=tol); fac(pos)=1; fac=1./fac;

dT_old = -((x.*fac)*(0:jmax)).*T_old;
dT_old(:,2:end) = dT_old(:,2:end) + (fac*(1:jmax)).*T_old(:,1:end-1);
dT_old(pos,:) = ones(length(pos),1)*(0:jmax).^2;

d2T_old = -fac*jall.^2.*T_old+(x(:,1).*fac)*ones(1,length(jall)).*dT_old;
d2T_old(pos,:) = ones(length(pos),1)*(jall.^2.*(jall.^2-1)/3);

[T,dT,d2T]=ClenshawCheb(x,jmax);
%
% G and Q have no singularities, but R has a cancellation involving singular
% terms. Do we compute it as (T_j-rT_j^\prime)/r^2?
%
Rdir = (T-x.*dT)./x2;
Rdir = Rdir(:,2:2:end);

jvec = 1:2:jmax;
jhalf = (jvec+1)/2;

nterm = jvec(end);
Rnew = (nterm-1)/factorial(nterm); % We need to add the correct sign here
for k=nterm-2:-2:1
    Rnew = -Rnew.*(x.^2*(jvec.^2-k^2)) + (k-1)/factorial(k);  
end
Rnew = Rnew.*(1./x*(jvec.*(-1).^jhalf));

if (ttype=='0')
    figure
    H=loglog(x,abs(Rnew),'-',x,abs(Rdir),'--');
    set(H,'LineWidth',1.5)
    for k=1:length(H)/2
        col=get(H(k),'Color');
        set(H(k+length(H)/2),'Color',col);
    end
    adjustFigProp
    xlabel('$r$','interpreter','latex')
    title('$(T_j(r)-rT_j^\prime(r))/r^2$','interpreter','latex')

    figure
    H = loglog(x,abs(Rdir-Rnew),'.');
    set(H,'MarkerSize',6)
    adjustFigProp
    xlabel('$r$','interpreter','latex')
    title('$(T_j(r)-rT_j^\prime(r))/r^2$','interpreter','latex')

    
elseif (ttype == '1')
    figure
    H = loglog(1-x,abs(dT-dT_old),'.');
    set(H,'LineWidth',1.5,'MarkerSize',6)
    adjustFigProp
    xlabel('$1-r$','interpreter','latex')
    title('Differences in $T_j^\prime(r)$','interpreter','latex')

    figure
    H = loglog(1-x,abs(d2T-d2T_old),'.');
    set(H,'LineWidth',1.5,'MarkerSize',6)
    adjustFigProp
    xlabel('$1-r$','interpreter','latex')
    title('Difference in $T_j^{\prime\prime}(r)$','interpreter','latex')
end

function adjustFigProp()
set(gca,'LineWidth',1.5,'FontSize',18)
set(gcf,'Color','w')
axis tight
grid
set(gca,'XTick',[1e-15 1e-10 1e-5 1e0])
set(gca,'YTick',[1e-15 1e-10 1e-5 1e0])
set(gca,'XMinorGrid','off')
set(gca,'YMinorGrid','off')
end

