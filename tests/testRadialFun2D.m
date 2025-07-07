close all
%

% We test the radial functions needed for 2D in RBF-QR here
%
jmax = 50;
n = 500;
x = logspace(-16,0,n)'; % Positive half line
x2 = x.^2;

[T,dT,d2T]=ClenshawCheb(x,jmax);
%
% As long as m > 0 F(r) in total has a positive power of r.
% if m = 0, p = 1 then we have T_j/r. This can only happen for odd j.
%
F = T./x;
F = F(:,2:2:end);
jvec = 1:size(F,2);
Fdiff = F(1,:)-(jvec*2-1).*(-1).^(jvec+1);
% Computes stably to very small values, but accuracy stops at 1e-14 approx.

figure
loglog(x,F)
title('F')
%
% G and Q have no singularities, but R has a cancellation involving singular
% terms. Do we compute it as (T_j-rT_j^\prime)/r^2?
%
R = (T-x.*dT)./x2;
R = R(:,2:2:end);
R(1,:)

jvec = 1:2:jmax;
jhalf = (jvec+1)/2;
R2 = (-1).^(jhalf+1)*1/3.*(jvec-1).*jvec.*(jvec+1);
max2 = max(abs(R2))
R2 = x*R2;
R2(:,1) = 0;

R3 = -(-1).^(jhalf+1)*1/30.*(jvec-3).*(jvec-1).*jvec.*(jvec+1).*(jvec+3);
max3 = max(abs(R3))
R3 = x.^3*R3;
R3(:,1:2) = 0;

R4 = (-1).^(jhalf+1)*1/840.*(jvec-5).*(jvec-3).*(jvec-1).*jvec.*(jvec+1).*(jvec+3).*(jvec+5);
max4 = max(abs(R4))
R4 = x.^5*R4;
R4(:,1:3)=0;

R5 = -(-1).^(jhalf+1)*1/45360.*(jvec-7).*(jvec-5).*(jvec-3).*(jvec-1).*jvec.*(jvec+1).*(jvec+3).*(jvec+5).*(jvec+7);
max5 = max(abs(R5))
R5 = x.^7*R5;
R5(:,1:4)=0;

j0 = jvec-1; j0(1)=1;
cpart = (jvec-1)./cumprod(jvec.*j0);
jpart = jmax*cumprod(jmax.^2-jvec.^2)
coeffmax = cpart(2:end).*jpart(1:end-1);
% This is validated that it reflects what is above.

%
% Now we want to find the breakoff point for the expansion.
% Where is the error of the expansion = tol?
% coeffmax*x.^jvec = tol
% x^j=tol/coeffmax
%
tol=1e-14;
xval = (tol./coeffmax).^(1./jvec(1:end-1))
%
% Then we also need the rounding error size.
%
xe = logspace(-16,-12,100)';
[Te,dTe,d2Te]=ClenshawCheb(xe,jmax);
Re = (Te-xe.*dTe)./xe.^2;
Re = Re(:,2:2:end);

for k=1:size(Re,2)
    pos = find(abs(Re(:,k))>0);
    Pcoeff(k,:) = polyfit(log10(xe(pos)),log10(abs(Re(pos,k))),1);
end
% Results are quite stable for different j
%
% a*x+b = log10(tol) => x = (log10(tol)-b)/a
%
xround = 10.^((log10(tol)-Pcoeff(:,2))./Pcoeff(:,1));
% Generally jmax also determines this one.

% Now xround says from which value direct computation is possible, while
% xval says how many components we need to reach that point.
% Say that The coefficient is always -1, then compute b from that.
%
% -log10(x) + b = log10(Re)
% A*b = log10(Re)+log10(x) mean(hl)
%
for k=1:size(Re,2)
    pos = find(abs(Re(:,k))>0);
    bcoeff(k) = mean(log10(xe(pos))+log10(abs(Re(pos,k))));
end
xround2 = 10.^(-(log10(tol)-bcoeff))
%Yes, in this way, the largest j determines the break point.
% This one can easily be computed as well.


% With the 5 terms, we can go up to 1e-12, choosing x=0.004 as the limit
% This would be for j up to 30. Let's now increase a bit to j=50. The error will be a bit larger or maybe much larger becasue of growth. We should stay away from where the error goes up. Probably we can see it from the size of the terms, where we should stay. Is the maximum (jmax^2-49)*(jmax^2-25)*(jmax^2-9)*(jmax^2-1)*jmax*max(x)^7

nterm = jvec(end);
R5b = (nterm-1)/factorial(nterm); % We need to add the correct sign here
for k=nterm-2:-2:1
    R5b = -R5b.*(x.^2*(jvec.^2-k^2)) + (k-1)/factorial(k);  
end
R5b = R5b.*(1./x*(jvec.*(-1).^jhalf));
% This is now correct apart from signs.Something seems to happen at x\approx 0.07, perhaps too few terms. At 0.05 1e-13 with terms up to 31
% We can precompute some terms x^2j^2 is one. x^2 in a matrix is another.

figure
loglog(x,abs(R),'-',x,abs(R2),'--')%,x,abs(R0(:,2:2:end)./x.^2),'o')
hold on
loglog(x,2*eps*jvec(:)*(1./x(:)'),'r-')
title('R')

figure
loglog(x,abs(R-R2))
title('Diff 2')

figure
loglog(x,abs(R3))
title('R3')

figure
loglog(x,abs(R-R2-R3))
title('Diff 3')

figure
loglog(x,abs(R-R2-R3-R4))
title('Diff 4')

figure
loglog(x,abs(R-R2-R3-R4-R5))
title('Diff 5')
%
% Stable to around 1e-8, the first is identically zero and should be computed as such. The limit should be zero, but all of them have the same slope.
% The limit does not look like zero. Can we compute what the difference should be with Clenshaw perhaps? We start with the case ep=0. T_j for j odd probably goes to zero in any case, and with ep, it will be different values.
%

P = ones(n,floor(jmax/2)+1);
for pp=1:size(P,2)-1
    P(:,pp+1) = x2.*P(:,pp);
end  
