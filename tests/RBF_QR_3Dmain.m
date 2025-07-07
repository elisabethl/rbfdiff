%--- Interpolates a function over the unit disc using RBF-QR and
%--- RBF-Direct. The best result in each case is displayed and the
%--- max-errors are displayed as functions of the shape parameter. 
clear
path(path,'/home/bette/Research/BFexp/MLAB') %%%%%%%%%%%%%LOCAL
N = 300;
Ne=10;
d=3;
%--- Approximately N Halton points in the unit sphere.
% Volume of a cube -1 to 1 is 2^3. Volume of sphere is pi r^3=pi.
xk=2*(halton(round(N/pi*8),d)-0.5);
%
[xk_p(:,3),xk_p(:,2),xk_p(:,1)]=cart2sph(xk(:,1),xk(:,2),xk(:,3));
xk_p(:,2) = pi/2-xk_p(:,2);

pos=find(xk_p(:,1)<=1); % Radius <= 1 in spherical coordinates
xk=xk(pos,:);  xk_p=xk_p(pos,:);

%--- Evaluation points from cartesian grid
t=-1:2/(Ne-1):1;
[x,y,z]=ndgrid(t,t,t);
x=x(:); y=y(:); z=z(:); 
pos = find(x.^2+y.^2+z.^2<=1);
x=x(pos);y=y(pos);z=z(pos);
[phi,th,r]=cart2sph(x,y,z);
th=pi/2-th;
xe_p=[r,th,phi];
xe = [x y z];

%--- Right hand side and exakt solution  
fun = @(x,y,z) sin(x.^2+2*y.^2-z.^2)-sin(2*x.^2+(y-0.5).^2+(z+0.1).^2);
f = fun(xk(:,1),xk(:,2),xk(:,3));  
u = fun(x,y,z);  
%--- Interpolate and evaluate
epvec=logspace(-2,log10(2),25);
[cx,xx]=meshgrid(xk(:,1));
[cy,yy]=meshgrid(xk(:,2));
[cz,zz]=meshgrid(xk(:,3));
r2 = (xx-cx).^2 + (yy-cy).^2 + (zz-cz).^2;
[cx,xx]=meshgrid(xk(:,1),x);
[cy,yy]=meshgrid(xk(:,2),y);
[cz,zz]=meshgrid(xk(:,3),z);
r2e = (xx-cx).^2 + (yy-cy).^2 + (zz-cz).^2;
for k=1:length(epvec)
  ep=epvec(k)
  E = RBF_QR_diffmat_3D(0,xk,ep,xe);  
  rbfqr(:,k) = E*f; %RBF_QR_3D(ep,xk_p,xe_p,f);
  A = exp(-ep^2*r2);
  B = exp(-ep^2*r2e);
  rbfdir(:,k) = B*(A\f);
  err_rbfqr(k)  = max(abs(rbfqr(:,k) -u));
  err_rbfdir(k) = max(abs(rbfdir(:,k)-u)); 
end
%--- Display solutions and errors
figure(1),clf
loglog(epvec,err_rbfqr,'b',epvec,err_rbfdir,'r--')
axis tight
xlabel('\epsilon')
ylabel('Maximum error')
legend('RBF-QR','RBF-Direct')

return
figure(2),clf
rw = [size(rr,1) 1:size(rr,1)];
xx=rr.*cos(tt); yy=rr.*sin(tt);
[mi,pos]=min(err_rbfqr);
solqr = reshape(rbfqr(:,pos),size(rr));
errqr = reshape(rbfqr(:,pos)-u,size(rr));
subplot(1,2,1)
H(1)=surf(xx(rw,:),yy(rw,:),solqr(rw,:)); axis square
title('RBF-QR: Best solution')
subplot(1,2,2)
H(2)=surf(xx(rw,:),yy(rw,:),errqr(rw,:)); axis square
title('RBF-QR: Smallest error')
set(H,'FaceColor','interp','EdgeColor','none')
figure(3),clf
[mi,pos]=min(err_rbfdir);
soldir = reshape(rbfdir(:,pos),size(rr));
errqr = reshape(rbfdir(:,pos)-u,size(rr));
subplot(1,2,1)
H(3)=surf(xx(rw,:),yy(rw,:),soldir(rw,:)); axis square
title('RBF-Direct: Best solution')
subplot(1,2,2)
H(4)=surf(xx(rw,:),yy(rw,:),errqr(rw,:)); axis square
title('RBF-Direct: Smallest error')
set(H,'FaceColor','interp','EdgeColor','none')

