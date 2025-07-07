function fe = Lpow_eval(T,re,ep,j,m,kvals)
% fe = Lpow_eval(T,re,ep,j,m,kvals)
% Evaluate Chebyshev series for the kth power of the Laplacian
% IN:
% T - coefficient matrix from Lpow_coef
% re - vector of evaluation points
% ep - shape parameter
% j - polynomial order
% m - expansion parameter, m = 0,...,floor(j/2)
% kvals - values for which to evaluate the kth Laplacian
% OUT:
% fe(*,:) - Delta^k f(r), where k are the values in kvals

re = re(:);
fe = zeros(length(re),length(kvals));
xmat = T{j+1}(:,:,m+1);

if max(kvals)>(size(T{1},1)-1)/2
    disp('Coefficients not computed, please run Lpow_coef again.');
    return
end

for i=1:length(kvals)
    K = kvals(i);
    % Start summing from high powers of epsilon
    for n=2*K+1:-1:1
        x = xmat(:,n+K^2);
        Pre = clenshaw(x,re);
        fe(:,i) = fe(:,i) + ep^(2*(n-1))*Pre;
    end
    fe(:,i) = fe(:,i).*exp(-ep^2*re.^2);
end

function s = clenshaw(A,x)
% Clenshaw's algorithm for evaluating Chebyshev series
N = length(A)-1;
n = length(x);

% U = zeros(n,3);
% idx = 1:3;
% for k=N:-1:1
%     U(:,idx(3)) = 2*x.*U(:,idx(2)) - U(:,idx(1)) + A(k+1);
%     idx = mod(idx,3) + 1;
% end
% s = A(1) - U(:,idx(1)) + x.*U(:,idx(2));

U = zeros(n,N+1);
U(:,N) = A(N+1)*ones(n,1);
for k = N:-1:2
    U(:,k-1) = 2*x.*U(:,k) - U(:,k+1) + A(k);
end
s = A(1) - U(:,2) + x.*U(:,1);

% U = zeros(n,2);
% U(:,1) = A(N+1)*ones(n,1);
% idx = [2 1];
% for k = N:-1:2
%     U(:,idx(1)) = 2*x.*U(:,idx(2)) - U(:,idx(1)) + A(k); 
%     idx = mod(idx,2) + 1;
% end
% s = A(1) - U(:,idx(1)) + x.*U(:,idx(2));