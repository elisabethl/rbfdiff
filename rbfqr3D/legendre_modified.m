function [y,sdy,y_over_s,sd2y,A,B] = legendre_modified(n,x,normalize)
%LEGENDRE Associated Legendre function.
%   P = LEGENDRE(N,X) computes the associated Legendre functions 
%   of degree N and order M = 0, 1, ..., N, evaluated for each element
%   of X.  N must be a scalar integer and X must contain real values
%   between -1 <= X <= 1.  
%
%   If X is a vector, P is an (N+1)-by-L matrix, where L = length(X).
%   The P(M+1,i) entry corresponds to the associated Legendre function 
%   of degree N and order M evaluated at X(i). 
%
%   In general, the returned array has one more dimension than X.
%   Each element P(M+1,i,j,k,...) contains the associated Legendre
%   function of degree N and order M evaluated at X(i,j,k,...).
%
%   There are three possible normalizations, LEGENDRE(N,X,normalize)
%   where normalize is 'unnorm','sch' or 'norm'.
%
%   The default, unnormalized associated Legendre functions are:
% 
%       P(N,M;X) = (-1)^M * (1-X^2)^(M/2) * (d/dX)^M { P(N,X) },
%
%   where P(N,X) is the Legendre polynomial of degree N. Note that
%   the first row of P is the Legendre polynomial evaluated at X 
%   (the M == 0 case).
%
%   SP = LEGENDRE(N,X,'sch') computes the Schmidt semi-normalized
%   associated Legendre functions SP(N,M;X). These functions are 
%   related to the unnormalized associated Legendre functions 
%   P(N,M;X) by:
%               
%   SP(N,M;X) = P(N,X), M = 0
%             = (-1)^M * sqrt(2*(N-M)!/(N+M)!) * P(N,M;X), M > 0
%
%   NP = LEGENDRE(N,X,'norm') computes the fully-normalized
%   associated Legendre functions NP(N,M;X). These functions are 
%   normalized such that
%
%            /1
%           |
%           | [NP(N,M;X)]^2 dX = 1    ,
%           |
%           /-1
%
%   and are related to the unnormalized associated Legendre
%   functions P(N,M;X) by:
%               
%   NP(N,M;X) = (-1)^M * sqrt((N+1/2)*(N-M)!/(N+M)!) * P(N,M;X)
%
%   Examples: 
%     1. legendre(2, 0.0:0.1:0.2) returns the matrix:
% 
%              |    x = 0           x = 0.1         x = 0.2
%        ------|---------------------------------------------
%        m = 0 |   -0.5000         -0.4850         -0.4400
%        m = 1 |         0         -0.2985         -0.5879
%        m = 2 |    3.0000          2.9700          2.8800  
%
%     2. X = rand(2,4,5); N = 2; 
%        P = legendre(N,X); 
%
%     so that size(P) is 3x2x4x5 and 
%     P(:,1,2,3) is the same as legendre(N,X(1,2,3)). 
%
%   Class support for input X:
%      float: double, single


%   Acknowledgment:
%
%   This program is based on a Fortran program by Robert L. Parker,
%   Scripps Institution of Oceanography, Institute for Geophysics and 
%   Planetary Physics, UCSD. February 1993.
%
%   Copyright 1984-2004 The MathWorks, Inc. 
%   $Revision: 5.22.4.4 $  $Date: 2004/07/05 17:02:05 $
%
%   Reference:
%     [1] M. Abramowitz and I.A. Stegun, "Handbook of Mathematical
%         Functions", Dover Publications, 1965, Ch. 8.
%     [2] J. A. Jacobs, "Geomagnetism", Academic Press, 1987, Ch.4.
%
%   Note on Algorithm:
%
%   LEGENDRE uses a three-term backward recursion relationship in M.
%   This recursion is on a version of the Schmidt semi-normalized 
%   Associated Legendre functions SPc(n,m;x), which are complex 
%   spherical harmonics. These functions are related to the standard 
%   Abramowitz & Stegun functions P(n,m;x) by
%
%       P(n,m;x) = sqrt((n+m)!/(n-m)!) * SPc(n,m;x)
%   
%   They are related to the Schmidt form given previously by:
%
%       SP(n,m;x) = SPc(n,0;x), m = 0
%                 = (-1)^m * sqrt(2) * SPc(n,m;x), m > 0
%   Also returns sdy which is needed in the partial derivative of Y
%   with respect to phi. This routine is called by RBFQR.
if nargin < 2
    error('Not enough input arguments')
elseif nargin > 5
    error('Too many input arguments')    
end

if numel(n) > 1 || ~isreal(n) || n < 0 || n ~= round(n)
    error('N must be a positive scalar integer');
end

if ~isreal(x) | max(abs(x)) > 1 | ischar(x)
    error('X must be real and in the range (-1,1)')
end

classin = superiorfloat(x);

% The n = 0 case
if n == 0 
    y = ones(size(x),classin);
    sdy = zeros(size(x),classin);
    sd2y = zeros(size(x),classin);
    y_over_s = zeros(size(x),classin); %Not true but it ends up being 0
    if nargin > 2 && strcmp(normalize,'norm')
        y = y/sqrt(2);
    end
    return
end

% Convert x to a single row vector
sizex = size(x);
x = x(:)';

rootn = sqrt(0:2*n);
s = sqrt(1-x.^2);
P = zeros(n+3,length(x),classin);
Q = zeros(n+3,length(x),classin);
R = zeros(n+3,length(x),classin);


% Calculate TWOCOT, separating out the x = -1,+1 cases first
twocot = x;

% Evaluate x = +/-1 first to avoid error messages for division by zero
k = find(x==-1);
twocot(k) = Inf(classin);

k = find(x==1);
twocot(k) = -Inf(classin);

k = find(s);
twocot(k) = -2*x(k)./s(k);

% Find values of x,s for which there will be underflow

sn = (-s).^n;
sn_minus_1 = -(-s).^(n-1);
sn_minus_2 = (-s).^(n-2);
tol = sqrt(realmin(classin));
ind = find(s>0 & abs(sn)<=tol);
if ~isempty(ind)
    % Approx solution of x*ln(x) = y 
    v = 9.2-log(tol)./(n*s(ind));
    w = 1./log(v);
    m1 = 1+n*s(ind).*v.*w.*(1.0058+ w.*(3.819 - w*12.173));
    m1 = min(n, floor(m1));

    % Column-by-column recursion
    for k = 1:length(m1)
        mm1 = m1(k);
        col = ind(k);
        P(mm1:n+1,col) = zeros(size(mm1:n+1))';

        % Start recursion with proper sign
        tstart = eps(classin);
        P(mm1,col) = sign(rem(mm1,2)-0.5)*tstart;
        if x(col) < 0
            P(mm1,col) = sign(rem(n+1,2)-0.5)*tstart;
        end

        % Recur from m1 to m = 0, accumulating normalizing factor.
        sumsq = tol;
        for m = mm1-2:-1:0
            P(m+1,col) = ((m+1)*twocot(col)*P(m+2,col)- ...
                  rootn(n+m+3)*rootn(n-m)*P(m+3,col))/ ...
                  (rootn(n+m+2)*rootn(n-m+1));
            sumsq = P(m+1,col)^2 + sumsq;
        end
        scale = 1/sqrt(2*sumsq - P(1,col)^2);
        P(1:mm1+1,col) = scale*P(1:mm1+1,col);
%         Q(1:mm1+1,col) = P(1:mm1+1,col)/s;
        Q(1:mm1+1,col) = P(1:mm1+1,col)./s(1:mm1+1).';
        R(1:mm1+1,col) = P(1:mm1+1,col)./(s(1:mm1+1).^2).';
        max(Q)
        max(R)
        pause
    end     % FOR loop
end % small sine IF loop

% Find the values of x,s for which there is no underflow, and for
% which twocot is not infinite (x~=1).

nind = find(x~=1 & abs(sn)>=tol);
if ~isempty(nind)

    % Produce normalization constant for the m = n function
    d = 2:2:2*n;
    c = prod(1-1./d);

    % Use sn = (-s).^n (written above) to write the m = n function
    P(n+1,nind) = sqrt(c)*sn(nind);
    P(n,nind) = P(n+1,nind).*twocot(nind)*n./rootn(end);
    
    Q(n+1,nind) = sqrt(c)*sn_minus_1(nind);
    Q(n,nind) = Q(n+1,nind).*twocot(nind)*n./rootn(end);
    
    R(n+1,nind) = sqrt(c)*sn_minus_2(nind);
    R(n,nind) = R(n+1,nind).*twocot(nind)*n./rootn(end);


    % Recur downwards to m = 0
    for m = n-2:-1:0
        P(m+1,nind) = (P(m+2,nind).*twocot(nind)*(m+1) ...
            -P(m+3,nind)*rootn(n+m+3)*rootn(n-m))/ ...
            (rootn(n+m+2)*rootn(n-m+1));
                
        Q(m+1,nind) = (Q(m+2,nind).*twocot(nind)*(m+1) ...
            -Q(m+3,nind)*rootn(n+m+3)*rootn(n-m))/ ...
            (rootn(n+m+2)*rootn(n-m+1));
        
        R(m+1,nind) = (R(m+2,nind).*twocot(nind)*(m+1) ...
            -R(m+3,nind)*rootn(n+m+3)*rootn(n-m))/ ...
            (rootn(n+m+2)*rootn(n-m+1)); 
    end
end % IF loop
    %
    % The first two should always be zero, no matter what $\mu$ is.
    %
    R(1:2,:) = 0;

y = P(1:n+1,:);
sdy= zeros(size(y));
sd2y= zeros(size(y));
A = zeros(size(y));
B = zeros(size(y));
y_over_s = Q(1:n+1,:);

% Polar argument   (x = +-1)
s0x1 = find(x == 1);
y(1,s0x1) = x(s0x1).^n;
sdy(2,s0x1) = 1/2*n*(n+1);
%
% General. There was some problems to find this bc, but now it shoud work
% for all n.
% 
R(3,s0x1) = 1/8*sqrt(prod(n-1:n+2));

s0xminus1 = find(x == -1);
y(1,s0xminus1) = x(s0xminus1).^n;
sdy(2,s0xminus1) = (-1)^(n)*1/2*n*(n+1);
%
% The same form for the other side but with some sign changes.
%
R(3,s0xminus1) = (-1)^n*1/8*sqrt(prod(n-1:n+2));

% y_over_s(2,s0x1)=1/2*n*(n+1);
y_over_s(2,s0x1)=-sqrt(n*(n+1))/2;
% y_over_s(2,s0xminus1)=(-1)^n*1/2*n*(n+1);
y_over_s(2,s0xminus1)=(-1)^n*sqrt(n*(n+1))/2;

if nargin < 3 || strcmp(normalize,'unnorm')
    % Calculate the standard A&S functions (i.e., unnormalized) by
    % multiplying each row by: sqrt((n+m)!/(n-m)!) = sqrt(prod(n-m+1:n+m))
    for m = 1:n-1
        y(m+1,:) = prod(rootn(n-m+2:n+m+1))*y(m+1,:);
        y_over_s(m+1,:) = prod(rootn(n-m+2:n+m+1))*y_over_s(m+1,:);
    end
    % Last row (m = n) must be done separately to handle 0!.
    % NOTE: the coefficient for (m = n) overflows for n>150.
    y(n+1,:) = prod(rootn(2:end))*y(n+1,:);
    y_over_s(n+1,:) = prod(rootn(2:end))*y_over_s(n+1,:);
    
%     nindex = find(x~=1 & x~=-1);
    nindex = 1:length(x);
    if ~isempty(nindex)
    
        sdy(n+1,nindex) = n*y(n,nindex);
        sdy(1,nindex) = -y(2,nindex);
        
        % Recur downwards to m = 0
        for m = 1:n-1
            sdy(m+1,nindex) = ((n+m)*(n-m+1)*y(m,nindex)-y(m+2,nindex))/2;
        end
    end % IF loop
elseif strcmp(normalize,'sch')
    % Calculate the standard Schmidt semi-normalized functions.
    % For m = 1,...,n, multiply by (-1)^m*sqrt(2)
    row1 = y(1,:);
    y = sqrt(2)*y;
    y(1,:) = row1;    % restore first row
    
    row1 = y_over_s(1,:);
    y_over_s = sqrt(2)*y_over_s;
    y_over_s(1,:) = row1;    % restore first row
    
    const1 = 1;
    for r = 2:n+1
        const1 = -const1;
        y(r,:) = const1*y(r,:);
        
        y_over_s(r,:) = const1*y_over_s(r,:);        
    end
elseif strcmp(normalize,'norm')
    % Calculate the fully-normalized functions.
    % For m = 0,...,n, multiply by (-1)^m*sqrt(n+1/2)
    y = sqrt(n+1/2)*y;
    
    y_over_s = sqrt(n+1/2)*y_over_s;
    
    R = sqrt(n+1/2)*R;
    
    const1 = -1;
    for r = 1:n+1
        const1 = -const1;
        y(r,:) = const1*y(r,:);
        
        y_over_s(r,:) = const1*y_over_s(r,:);
        
        R(r,:) = const1*R(r,:);
    end
%     nindex = find(x~=1 & x~=-1);

%
% Note that all points are included here. That is, it does not matter if
% we set the end points correctly before, since they are not used again.
%
    nindex = 1:length(x);

    if ~isempty(nindex)
        sdy(n+1,nindex) = -sqrt(n/2)*y(n,nindex);
        sdy(1,nindex) = sqrt(n*(n+1))*y(2,nindex);
%         sdy(1,nindex) = n*(n+1)*y(2,nindex);

        % Recur downwards to m = 0
        for m = 1:n-1
            sdy(m+1,nindex) = -(sqrt((n+m)*(n-m+1))*y(m,nindex)-...
                sqrt((n+m+1)*(n-m))*y(m+2,nindex))/2;
        end
        
        sd2y(n+1,nindex) = -sqrt(n/2)*sdy(n,nindex);
        sd2y(1,nindex) = sqrt(n*(n+1))*sdy(2,nindex);
        for m = 1:n-1
            sd2y(m+1,nindex) = -(sqrt((n+m)*(n-m+1))*sdy(m,nindex)-...
                sqrt((n+m+1)*(n-m))*sdy(m+2,nindex))/2;
        end
        %
        % We compute the combinations that are needed for the second
        % derivatives of the spherical harmonics in Y.m
        %
        if (n==1)
          A(n+1,nindex) = -y(2,:);
          B(n+1,nindex) = -y(2,:);
        else
          A(n+1,nindex) = -n^2*R(n+1,nindex) ...
              +x(nindex)*sqrt(0.5*n).*y_over_s(n,nindex);
          B(n+1,nindex) = -n*(R(n+1,nindex) ...
                          -x(nindex).*sqrt(0.5*n).*y_over_s(n,nindex));
        end
        A(1,nindex) = -x(nindex).*sqrt(n*(n+1)).*y_over_s(2,nindex);
        %
        B(1,nindex) = 0;
        for m=1:n-1
          A(m+1,nindex) = -m*(m-1)*R(m+1,nindex)-m*y(m+1,nindex) ...
              -rootn(n+m+1+1)*rootn(n-m+1)*x(nindex).*y_over_s(m+2,nindex);
          %
          B(m+1,nindex) = m*(m-1)*R(m+1,nindex)-m*m*y(m+1,nindex) ...
              -m*rootn(n+m+1+1)*rootn(n-m+1)*x(nindex).*y_over_s(m+2,nindex);
        end  
    end % IF loop
else
    error(['Normalization option ''' normalize ''' not recognized'])
end

% y_over_s(2,s0x1)=-1/2*n*(n+1);
% y_over_s(2,s0x1)=1/2*n*(n+1);
% y_over_s(2,s0xminus1)=(-1)^n*1/2*n*(n+1); % sign here?

% Restore original dimensions.
if length(sizex) > 2 || min(sizex) > 1
    y = reshape(y,[n+1 sizex]);
    sdy = reshape(sdy,[n+1 sizex]);
    y_over_s = reshape(y_over_s,[n+1 sizex]);
end