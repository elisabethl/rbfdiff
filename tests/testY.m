%
% Test the derivatives of all spherical harmonics Y_mu^nu for a fixed mu
%
function testY(mu,N)
%
% First test the derivatives wrt to theta for a fixed value of fi
%
fi = 2*pi*rand(1);
%
% We use N values to evaluate the derivatives numerically using finite
% differences and then compare with the computed values at the midpoints.
%
theta = linspace(0,pi,N)';
h = (theta(2)-theta(1));
theta2 = theta(2:N)-h/2; % The shifted values
npi = 1;
%
% We measure the error away from the nearly singular expressions, where the
% values computed for comparison are not accurate.
%
tol = 1e-3; 
pos = find(abs(theta)>tol & abs(pi-theta)>tol);
pos2 = find(abs(theta2)>tol & abs(pi-theta2)>tol);

[YY,Yth,Yfi_over_s,Yfi,Ythth,Ythfi,Yfifi]=Y(mu,theta,fi);
[YY2,Yth2,Yfi_over_s2,Yfi2,Ythth2,Ythfi2,Yfifi2]=Y(mu,theta2,fi);

dim = 1;
%
% The derivative wrt to theta
%
numder = (1/h)*diff(YY,1,dim);

figure(1),clf
numder = reshape(numder,N-1,(mu+1)^2);
Yth = reshape(Yth,N,(mu+1)^2);
H = plot(theta,Yth,theta(2:N)-h/2,numder,'.');
xlabel('$\theta$','interpreter','latex');
fixColor(H)
adjustFigProp(npi);
title('$\partial Y/\partial\theta$','interpreter','latex')

Yth2 = reshape(Yth2,N-1,(mu+1)^2);
for k=1:(mu+1)^2
  norm_th(k) = norm(Yth2(pos2,k)-numder(pos2,k))/sqrt(length(pos2));
end  
norm_th
%
% The second derivative wrt to theta
%
numder = (1/h)*diff(Yth,1,dim);
figure(2),clf
numder = reshape(numder,N-1,(mu+1)^2);
Ythth = reshape(Ythth,N,(mu+1)^2);
H=plot(theta,Ythth,theta2,numder,'.');
fixColor(H)
adjustFigProp(npi);
xlabel('$\theta$','interpreter','latex')
title('$\partial^2 Y/\partial\theta^2$','interpreter','latex')

for k=1:(mu+1)^2
  norm_thth(k) = norm(Ythth2(pos2,k)-numder(pos2,k))/sqrt(length(pos2));
end  
norm_thth
%
% Then test the derivative wrt to fi for a fixed value of theta
%
theta = pi/3*ones(N,1); % Choosing a safe value
%
% We use N values to evaluate the derivatives numerically and then compare
%
fi = linspace(0,2*pi,N)';
h = (fi(2)-fi(1));
npi = 2;
[YY,Yth,Yfi_over_s,Yfi,Ythth,Ythfi,Yfifi]=Y(mu,theta,fi);
%
% We also need the values at the midpoints 
%
[YY2,Yth2,Yfi_over_s2,Yfi2,Ythth2,Ythfi2,Yfifi2]=Y(mu,theta(1:N-1),fi(2:end)-h/2);

isint = 1/sin(theta(1)); 
ct = cos(theta(1));

numder = (1/h)*diff(YY,1,1);

numfifi = (1/h)*diff(Yfi_over_s,1,dim);
numfifi= isint*(numfifi+ct*Yth2);

numthfi = (1/h)*diff(Yth,1,1);
numthfi = isint*(Yfi_over_s2-ct*numthfi);

figure(3),clf
numder_s = isint*reshape(numder,N-1,(mu+1)^2);
Yfi_s = reshape(Yfi_over_s,N,(mu+1)^2);
H = plot(fi,Yfi_s,fi(2:N)-h/2,numder_s,'.');
fixColor(H)
adjustFigProp(npi);
xlabel('$\varphi$','interpreter','latex')
title('$\sin^{-1}\theta\partial Y/\partial\varphi$','interpreter','latex')
%

figure(4),clf
numder = reshape(numder,N-1,(mu+1)^2);
Yfi = reshape(Yfi,N,(mu+1)^2);
H = plot(fi,Yfi,fi(2:N)-h/2,numder,'.');
fixColor(H)
adjustFigProp(npi);
xlabel('$\varphi$','interpreter','latex')
title('$\partial Y/\partial\varphi$','interpreter','latex')


figure(5),clf
Yfifi = reshape(Yfifi,N,(mu+1)^2);
numfifi = reshape(numfifi,N-1,(mu+1)^2);
H = plot(fi,Yfifi,fi(2:N)-h/2,numfifi,'.');
fixColor(H)
adjustFigProp(npi);
xlabel('$\varphi$','interpreter','latex')
title('$A$','interpreter','latex')

Yfifi2 = reshape(Yfifi2,N-1,(mu+1)^2);
diffA = reshape(max(abs(Yfifi2-numfifi)),mu+1,mu+1)

figure(6),clf
Ythfi = reshape(Ythfi,N,(mu+1)^2);
numthfi = reshape(numthfi,N-1,(mu+1)^2);
H = plot(fi,Ythfi,fi(2:N)-h/2,numthfi,'.');
fixColor(H)
adjustFigProp(npi);
xlabel('$\varphi$','interpreter','latex')
title('$B$','interpreter','latex')

Ythfi2 = reshape(Ythfi2,N-1,(mu+1)^2);
diffB = reshape(max(abs(Ythfi2-numthfi)),mu+1,mu+1)
%
%
function fixColor(H)
np = length(H);
for k=1:np/2
  col=get(H(k),'Color');
  set(H(k+np/2),'Color',col);
end 

function adjustFigProp(npi)
set(gca,'LineWidth',1.5,'FontSize',18)
set(gcf,'Color','w')
axis tight
grid
if (npi==1)
    set(gca,'XTick',[0 pi/2 pi],'TickLabelInterpreter','latex','XTickLabel',{0,'$\pi/2$','$\pi$'})
else
    set(gca,'XTick',[0 pi 2*pi],'TickLabelInterpreter','latex','XTickLabel',{0,'$\pi$','$2\pi$'})
end
