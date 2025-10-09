function [Psi,Rt]=InitPsi_3D(ep,xk)
% The purpose of this function is to compute \tilde{R} that is used
% for evaluating the basis functions \Psi used in the RBF-QR method.
%--- ep (scalar) : The shape parameter
%--- xk(1:N,1:3) : The center points in spherical coordinates and
%                  scaled to the unit ball. 
N = size(xk,1);
mp = eps; % machine precision
tol = 1e4*mp;
tolqr = 2;
%--- Find out how many columns are needed at the least to include N
jN = degree(N,3);

% Figure out what indices we need for the 3D-case. We order the
% functions as j first then m, then nu. I think there may be many
% zeros, but this is something we can perhaps gradually improve as
% it goes.
[Psi,C,T] = AddCBlocks([],[],[],ep,xk,jN);
M = length(Psi.j);
%--- QR-factorize the coefficient matrix and compute \tilde{R}
%--- incrementally block by block, until we find a block that is
%--- insignificant in absolute value compared to the previous ones.
[Q,R]=IncQR(C,ep^2,Psi.j,tolqr);

true = 1;
q = Q.q;
posN = -1;
jmax = jN;
%R = zeros(N,M); columns(1)=1; lindep=zeros(0,1); q=2;
%[Q,R(:,1)] = qr(C(:,1));
%for p=2:M
while true
    %--- Add one block to C
    jmax = jmax+1;
    [Psi,C,T,Cnew] = AddCBlocks(Psi,C,T,ep,xk,jmax);

    %--- If we have less than N columns, then we cannot be done
    if (q > N)
        if (posN<0) % First time, find out location to compare with
	    jN = Psi.j(R.order(N));
	    pos = find(Psi.j==jN);
	    posN = pos(1); % The location of d_j,0
        end
      
        %--- Check the magnitude of the latest block. 
        jtest = Psi.j(end);
        pos = find(Psi.j==jtest);
        pos = pos(1);
        p1 = [1 posN]; % First and Nth, either one may be smallest
        p2 = [pos];
        relsc = EvalD_3D(ep,p1,p2,Psi.j,Psi.m,Psi.p);
        relsc = max(relsc);

        %--- Block small enough means we are done.
        if (relsc*exp(0.223*jtest-0.012-0.649*mod(jtest,2)) < mp)
	    break
        end
    end  
    %--- Update the QR-factorization with the new block
    pos = find(Psi.j==Psi.j(end));
    [Q,R]=IncQR(Cnew,ep^2,Psi.j(pos),tolqr,Q,R);
    q = Q.q;
end
    

M = length(R.order);
Rt = R.R(1:N,1:N)\R.R(1:N,N+1:M);

p1 = R.order(1:N);    p2 = R.order((N+1):M);
 
if (M>N)
    D = EvalD_3D(ep,p1,p2,Psi.j,Psi.m,Psi.p);
    Rt = D.*Rt;
end
Psi.ep = ep;
Psi.xk = xk;
Psi.Rt = Rt;
Psi.columns = R.order;


function [Psi,C,T,C1] = AddCBlocks(Psi,C,T,ep,xk,jmax)
%
% The first time, Psi should be empty, then it contains previous
% results. The same goes for C.
% jmax is optional, if not present, just one block is added.
% if jmax is present, all blocks up to jmax are added.
% This could be useful since we know initially the least number of
% blocks we need.
%
if isempty(Psi) % Determine starting point for new columns
    j0 = 0;
    T.Yk = [];
    T.Pk = ones(size(xk,1),jmax+1); % Fill with ones up to the first end
    T.rscale = exp(-ep^2*xk(:,1).^2); % Row scaling of C = exp(-ep^2*r_k^2)
    Psi.j = zeros(0,1);
    Psi.m = zeros(0,1);
    Psi.p = zeros(0,1);
    Psi.nu = zeros(0,1);
    C = zeros(size(xk,1),0);
else
    j0 = Psi.j(end) + 1;
end
if (nargin==5) % Determine final point
    jmax = j0;
end
j = zeros(0,1); m=zeros(0,1); p=zeros(0,1); odd=mod(j0+1,2);

nu = zeros(0,1);
for k=j0:jmax
    odd = 1-odd;
    nk = (k+1)*(k+2)/2;
    j = [j; k*ones(nk,1)];   p = [p; odd*ones(nk,1)];
    for mm=0:(k-odd)/2
        nn = -(2*mm+odd):(2*mm+odd);
        nu = [nu; nn'];
        m = [m; mm*ones(size(nn'))];
    end
end

T.Yk = Y(jmax,xk(:,2),xk(:,3),T.Yk);

for k=max(1,j0):jmax
    T.Pk(:,k+1) = xk(:,1).*T.Pk(:,k);
end
%--- Compute the coefficient matrix
M = length(j);

cscale = ones(1,M);
pos = find(nu == 0);
cscale(pos) = 0.5*cscale(pos);
pos = find(j-2*m == 0);
cscale(pos) = 0.5*cscale(pos);

C1 = T.Pk(:,j+1); % The powers of r_k and then the trig part

% This part is also inefficient, should be able to get rid of loop
% Not going to bother for now
for k=1:M
    mu = 2*m(k)+p(k);
    if (nu(k)>=0)
        C1(:,k) = C1(:,k).*T.Yk(:,mu-nu(k)+1,mu+1);
    else
        C1(:,k) = C1(:,k).*T.Yk(:,mu+1,mu+nu(k)+1);
    end
end

C1 = C1.*(T.rscale*cscale);

a = [(j-2*m+1)/2 (j-2*m+2)/2];
b = [j-2*m+1 (j-2*m-p+2)/2 (j+2*m+p+3)/2];

z = ep.^4*xk(:,1).^2;

for k=1:M
    C1(:,k) = C1(:,k).*hypergeom23(a(k,:),b(k,:),z);
end
Psi.j = [Psi.j; j];
Psi.m = [Psi.m; m];
Psi.p = [Psi.p; p];

Psi.nu = [Psi.nu; nu];
C = [C C1];



