%
% m is the order of the legendre polynomials
% N is the number of points used for the numerical derivatives
%
function testLegendre(m,N)
 
theta = linspace(0,pi,N);  
z = cos(theta);
[y,sdy,y_s,sd2y,A,B] = normalizedLegendre(m,z);
%
% Compute the first and second derivatives numerically
%
h = theta(2)-theta(1);
dim = 2;
yprime = diff(y,1,dim)/h;
ybis = diff(y,2,dim)/h^2;

figure(1),clf
H=plot(theta,sdy,theta(2:end)-h/2,-yprime,'.');
np = length(H)
for k=1:np/2
  col=get(H(k),'Color');
  set(H(k+np/2),'Color',col);
end  
adjustFigProp
title('$\frac{\partial P_\mu^\nu}{\partial\theta}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')


figure(2),clf
H=plot(theta,sd2y,theta(2:end-1),ybis,'.');
for k=1:np/2
  col=get(H(k),'Color');
  set(H(k+np/2),'Color',col);
end
adjustFigProp
title('$\frac{\partial^2 P_\mu^\nu}{\partial\theta^2}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')
%
% We next check if y_s=y/sin(theta)
% Note that the value for nu=0 is unboundad and is not used in the code
%
Q=y(2:end,2:end-1)./sin(theta(2:end-1));

figure(3),clf
H=plot(theta,y_s(2:end,:),theta(2:end-1),Q,'.');
for k=1:length(H)/2
  col=get(H(k),'Color');
  set(H(k+length(H)/2),'Color',col);
end
adjustFigProp
title('$\frac{P_\mu^\nu}{\sin\theta}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')

%
% A_mu^\nu = -nu^2y/sin^2theta-cos(theta)/sin(theta)*sdy
% 
ind = 2:length(theta)-1; % sin(theta) is zero at the end points
isin = 1./sin(theta(ind)); 
for nu=0:m
  Anum(nu+1,:) = -nu^2*y(nu+1,ind).*isin.^2-z(ind).*isin.*sdy(nu+1,ind);
end
figure(4),clf
H=plot(theta,A(1:end,:),theta(ind),Anum(1:end,:),'.');
for k=1:length(H)/2
  col=get(H(k),'Color');
  set(H(k+length(H)/2),'Color',col);
end
adjustFigProp
title('$A_\mu^\nu=-\frac{\nu^2}{\sin^2\theta}P_\mu^\nu-\frac{\cos\theta}{\sin\theta}\frac{\partial P_\mu^\nu}{\partial\theta}$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')

%
% B_mu^\nu = -nu( (1+cos^2 theta)y/sin^2theta+2*cos(theta)/sin(theta)*sdy
% 
for nu=0:m
  Bnum(nu+1,:) = -nu*(y(nu+1,ind).*isin.^2+z(ind).*isin.*sdy(nu+1,ind));
end
figure(5),clf
% The first one is identically zero
H=plot(theta,B(2:end,:),theta(ind),Bnum(2:end,:),'.');
for k=1:length(H)/2
  col=get(H(k),'Color');
  set(H(k+length(H)/2),'Color',col);
end
adjustFigProp
title('$B_\mu^\nu=-\nu\left(\frac{1}{\sin^2\theta}P_\mu^\nu+\frac{\cos\theta}{\sin\theta}\frac{\partial P_\mu^\nu}{\partial\theta}\right)$','interpreter','latex')
xlabel('$\theta$','interpreter','latex')

function adjustFigProp()
set(gca,'LineWidth',1.5,'FontSize',18)
set(gcf,'Color','w')
axis tight
grid
set(gca,'XTick',[0 pi/2 pi],'TickLabelInterpreter','latex','XTickLabel',{0,'$\pi/2$','$\pi$'})



