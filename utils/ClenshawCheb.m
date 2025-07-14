function [T,dT,d2T] = ClenshawCheb(x,jmax)
x = x(:); % We can only work with a vector of 1D-values
  n = size(x,1);
  v = zeros(n,jmax+1);
  v(:,1) = 1;
  v(:,2) = 2*x;
  T = zeros(n,jmax+1);
  T(:,1) = 1; %0.5;
  T(:,2) = x;
  for j=2:jmax
    v(:,j+1) = 2*x.*v(:,j) - v(:,j-1);
    T(:,j+1) =   x.*v(:,j) - v(:,j-1);
  end

  w = zeros(n,jmax+1);
  w(:,2) = 2*v(:,1);
  dT = zeros(n,jmax+1);
  dT(:,2) = 1;
  R = zeros(n,jmax+1);
  for j=2:jmax
     w(:,j+1) = 2*x.*w(:,j) - w(:,j-1) + 2*v(:,j);
     dT(:,j+1) =   x.*w(:,j) - w(:,j-1) +   v(:,j);
     %R(:,j+1) =  -v(:,j-1)+ x.*w(:,j-1) - x.^2.*w(:,j);
  end

  z = zeros(n,jmax+1);
  d2T = zeros(n,jmax+1);
  for j=2:jmax
        z(:,j+1) = 2*x.*z(:,j) - z(:,j-1) + 4*w(:,j);
      d2T(:,j+1) =   x.*z(:,j) - z(:,j-1) + 2*w(:,j);
  end
  % x*v(:,j)-v(:,j-1)-x^2*w(:,j)+x*w(:,j-1)-x*v(:,j) =
  % -v(:,j-1)-x^2*w(:,j)+x*w(:,j-1)
  % This is unstable from j=24 and upp for x below 1e-9. More stable than the
  % other formulation though.
  
