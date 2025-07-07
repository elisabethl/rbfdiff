clear
%--- The RBF-QR method is used for approximating a function and its 
%--- first and second derivative (in one dimension)
fnum = 1;

if (fnum==1)
    fun = @(x) ones(size(x));
    funx = @(x) zeros(size(x));
    funxx = @(x) zeros(size(x));
elseif (fnum==2)    
    fun = @(x) 25./(25+(x-0.2).^2);
    funx = @(x) -50*(x-0.2)./(25+(x-0.2).^2).^2;
    funxx = @(x) 200*(x-0.2).^2./(25+(x-0.2).^2).^3- 50./(25+(x-0.2).^2).^2;
elseif (fnum==3)    
    fun = @(x) exp(-(x-0.1).^2);
    funx = @(x) -2*(x-0.1).*exp(-(x-0.1).^2);
    funxx = @(x) -2*exp(-(x-0.1).^2) + 4*(x-0.1).^2.*exp(-(x-0.1).^2);
elseif (fnum==4)    
    fun = @(x) sin(x.^2)-sin(2*x.^2);
    funx = @(x) 2*x.*cos(x.^2)-4*x.*cos(2*x.^2);
    funxx = @(x) 2*cos(x.^2)-4*cos(2*x.^2)-4*x.^2.*sin(x.^2) + 16*x.^2.*sin(2*x.^2);
end

%--- Choose problem size and range of epsilon
N =  20; % The number of node points
Ne= 400;  % The number of evaluation points

epvec=logspace(-2,log10(4),50);

%--- Generate node points. Here, we choose a Chebyshev grid, which has
%--- been shown to be beneficial in the small epsilon range.
theta = (pi:-pi/(N-1):0)';
xk = cos(theta);

%--- Uniform evaluation points
xe = linspace(-1,1,Ne)';
xe = linspace(1-1e-6,1,Ne)';
%xe = linspace(1-1e-6,1,Ne)'; 
xe = linspace(1,1+1e-2,Ne)'; % Up to here we get decent results

%--- Right hand side and exact solution  
f = fun(xk);  
u = fun(xe); 
ux = funx(xe); 
uxx = funxx(xe); 

%--- Distance matrices for the RBF-Direct approach, which will be used
%--- for comparison
[cc,xx]=meshgrid(xk);
rk2 = (xx-cc).^2;
[cc,xx]=meshgrid(xk,xe);
re2 = (xx-cc).^2;
ndiff = 1.5; % Compute operators up to Laplacian
for k=1:length(epvec)
  ep=epvec(k);
  %--- Compute the interpolant using RBF-Direct
  A = exp(-ep^2*rk2);
  B = exp(-ep^2*re2);
  rbfdir(:,k) = B*(A\f);
  %--- Compute the interpolant and its derivatives using RBF-QR
  [E,Psi,Te] = RBF_QR_diffmat_1D(0,xk,ep,xe);
  G = RBF_QR_diffmat_1D(1,Psi,Te);
  L = RBF_QR_diffmat_1D(1.5,Psi,Te);
  rbfqr(:,k) = E*f;
  rbfx(:,k)  = G{1}*f;
  rbfxx(:,k) = L*f;
  %--- Use the other method of computation
  C = 0; R = 1;
  Psi = RBF_QR_diffmat_1D(ndiff,xk,ep,C,R);
  [E,Te] = RBF_QR_diffmat_1D(0,Psi,xe);
  G = RBF_QR_diffmat_1D(1,Psi,Te);
  L = RBF_QR_diffmat_1D(1.5,Psi,Te);
  rbfqr2(:,k) = E*f;
  rbfx2(:,k) = G{1}*f;
  rbfxx2(:,k) = L*f;
  %--- Compute the errors
  errdir(k) = max(abs(rbfdir(:,k)-u)); 
  errqr(k) = max(abs(rbfqr(:,k) -u));
  errx(k) = max(abs(rbfx(:,k)-ux)); 
  errxx(k) = max(abs(rbfxx(:,k)-uxx));
  errqr2(k) = max(abs(rbfqr2(:,k) -u));
  errx2(k) = max(abs(rbfx2(:,k)-ux)); 
  errxx2(k) = max(abs(rbfxx2(:,k)-uxx)); 
end
%--- Display the errors
figure(1),clf
loglog(epvec,errqr,'b',epvec,errdir,'r--',epvec,errx,'k',epvec,errxx,'m')
hold on
loglog(epvec,errqr2,'b--',epvec,errx2,'k--',epvec,errxx2,'m--')
axis tight
xlabel('\epsilon','FontSize',15)
ylabel('Maximum error','FontSize',15)
legend('RBF-QR','RBF-Direct','1st deriv','2nd deriv')

%--- Find the best solutions in both cases
figure(2),clf
[miqr,posqr]=min(errqr);

subplot(1,2,1)
H(1)=plot(xe,rbfqr(:,posqr)); 
title('RBF-QR: Best solution')
subplot(1,2,2)
H(2)=plot(xe,rbfqr(:,posqr)-u); 
title('RBF-QR: Smallest error')

figure(3),clf
[mid,posd]=min(errdir);

subplot(1,2,1)
H(3)=plot(xe,rbfdir(:,posd)); 
title('RBF-Direct: Best solution')
subplot(1,2,2)
H(4)=plot(xe,rbfdir(:,posd)-u); 
title('RBF-Direct: Smallest error')

figure(4),clf
plot(xe,rbfx(:,posqr),xe,funx(xe),'--')
title('Best gradient with exact')
legend('RBF-QR','Exact')

figure(5),clf
plot(xe,rbfxx(:,posqr),xe,funxx(xe),'--')
title('Best 2nd derivative with exact')
legend('RBF-QR','Exact')
