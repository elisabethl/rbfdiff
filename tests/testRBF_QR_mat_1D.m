%
% The first and second derivatives have removable singularities at 1 and -1.
% This test is to ensure that we get correct results for values close to and at
% these points.
%
dim = 1;
ep = 0.1;
%
% We use Chebyshev points to have a numerically stable approximation near the
% boundary
%
N = 30;
theta = (pi:-pi/(N-1):0)';
xk = cos(theta);
%
% We put the evaluation point close to -1
%
xdiff = logspace(-16,-3,200)';
xe = -1 + xdiff;

[Psi] = InitPsi_1D(ep,xk);
Dx = RBF_QR_mat_1D(Psi,'x',xe);
Dxx = RBF_QR_mat_1D(Psi,'xx',xe);

figure(1),clf
semilogy(xdiff,abs(Dx(:,15)))
% First derivative looks ok up to 5e-15. There is a slow variation that should be captured.

figure(2),clf
semilogy(xdiff,abs(Dxx(:,15)))
ax=axis;
ax(3:4) = 3e5*[-1 1];
axis(ax)
% Here errors are very large, seems to stabilize around 1e-8. We need to look at both the limit and the values in between.

%
% We may also want to check extrapolation after this deed is done.
%


%
% We may also want to check the derivative results for a function
%
