function [Psi]=InitPsi_1D(ep,xk)
% The purpose of this function is to compute \Psi, which is a
% representation of the RBF-QR basis functions. 
% Indata:  
%--- ep (scalar) : The shape parameter
%--- xk(1:N) : The center points

  N = size(xk,1); 
  mp = eps; % machine precision

%--- The sign is extracted and handled separately 
  rk = abs(xk);
  sk = ones(size(xk));
  pos = find(xk<0);
  sk(pos)=-1; % Can be interpreted as the 1-D angle
  
%--- Find out how many columns are needed. The hypergeometric function
%--- increases the value at j=0 compared with the rest, so this
%--- estimate should be overly pessimistic
  jN=N-1;
  jmax = 1; ratio=ep^2;
  % Find out if d_00 or d_jN is smallest and compute dj
  while (jmax<jN & ratio > 1) 
    jmax = jmax+1;
    ratio = ep^2/jmax*ratio;
  end
  
  if (ratio < 1) % d_jN was smallest
    jmax = jN; ratio=1;
  end	
  
  ratio = ep^2/(jmax+1); 
  while (ratio > mp) 
    jmax = jmax + 1;
    ratio = ep^2/(jmax+1)*ratio;
  end  
  jmax = max(jmax,jN); % Must always have at least N columns
  
%--- Compute the coefficient matrix C, up to degree jmax 
  [C,hs] = ComputeC(ep,rk,sk,jmax);

%--- QR-factorize the coefficient matrix and compute \tilde{R}
  [Q,R] = qr(C); 
  M = jmax+1;
  Rt = R(1:N,1:N)\R(1:N,N+1:M);
  p1 = (1:N);    p2 = (N+1):M;   
 
  if (M>N) % If there is a part to scale
    D = EvalD_1D(ep,p1,p2,(0:jmax)');
    Rt = D.*Rt;
    % Do the scaling from the hypergeometric 
    d1 = spdiags(hs(1:N),0,N,N); % Inverse of 1/hs
    d2 = spdiags(1./hs(N+1:M),0,M-N,M-N);
    Rt = d1*Rt*d2;
  end  
  
  Psi.ep = ep;
  Psi.xk = xk;
  Psi.Rt = Rt;
  Psi.M = M;
  Psi.jmax = jmax;

function [C,hs] = ComputeC(ep,rk,sk,jmax)
% Compute the coefficient matrix C up to degree jmax
%--- ep (scalar) : The shape parameter
%--- rk(1:N) : The absolute value of the center points 
%--- sk(1:N) : The sign of the center coordinates
%--- jmax (scalar) : The degree to stop at
 
  [ma,pos] = max(rk);  
  Pk = ones(size(rk,1),jmax+1); % Fill with ones up to the first end
  for k=1:jmax
      Pk(:,k+1) =rk(:,1).*Pk(:,k);
  end
  
%--- Compute the new blocks of the coefficient matrix
  rk2 = rk.*rk;
  M = jmax+1; 
  rscale = exp(-ep^2*rk2); % Row scaling of C = exp(-ep^2*r_k^2)
  cscale = 2*ones(1,M);    % Column scaling is just the constant 2
  C = Pk(:,1:M);           % The powers of r_k
  C = C.*(rscale*cscale);
  
  b = (0:jmax)'+1;
  z = ep.^4*rk2;  
  for k=1:M
    v = hypergeom01(b(k),z);
    hs(k,1)=1/v(pos);
    C(:,k) = C(:,k).*v*hs(k);
  end  

%-- For each odd j, multiply by the sign vector
  for k=2:2:M
    C(:,k) = C(:,k).*sk;
  end  



