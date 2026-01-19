function [phi]=bmp(rho,r,nprime,dim)
%
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% DIM(1:2) are the dimensions for the mixed second derivative if nprime
% is 'm2'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
% We assume that the r that comes in is unscaled and rho is the unscaled radius of the basis functions.
%
% Assume a one-dimensional problem if no dimension is given
% 
if nargin<=3
  dim=1;
end
%
% For the case of L or L2 operators, we need to know the number of dimensions
%
if (nprime(1)=='L')
  if (size(r,3)==1)
    nd=dim;
  else
    nd=size(r,3)-1;
  end
end
%
% For the mixed derivative, the dimensions must be given for all d
%
if (nprime(1)=='m')
  p = str2num(nprime(2));
  if (length(dim)~=p)
    error('For the mixed derivative mp, dim=dim(1:p)')
  elseif (all(dim==dim(1)))
    error('For mixed derivatives, dim(:) cannot have only one value')
  end  
end

phi=zeros(size(r,1),size(r,2));

% Here we scale the given radius with rho, such that we can compute all derivatives in the reference scaling. Later we add the scaling from differentiation.
r = (1/rho)*r;

% Make sure every value for r >= 1 becomes zero
r0 = sq(r(:,:,1));
rIsIn = (r0<1);
r0 = r0.*rIsIn; % Make r = zero from 1 and out
tmpPw = (1./(1-r0.^2)).*rIsIn; % This result also needs to be zero outside 
tmp =   exp(-tmpPw).*rIsIn;    % To put zero at r=1 and outside

% Precompute derivatives up to the order we need
switch nprime
  case {'0','1','2','3','4'}
    n = str2num(nprime);
  case {'L'}
    n = 2;
  case {'m2','m3','m4'}
    n = str2num(nprime(2));
  case {'L2'}
    n = 4;
  otherwise
    error('Input nprime not allowed to function bmp')
end

fac = tmp.*tmpPw.^2;
if (n >= 1)
  f1 = -fac;
end
if (n >= 2)
  fac = tmpPw.^2.*fac;
  f2 = (2*r0.^2-1).*fac;
end
if (n >= 3)
  fac = tmpPw.^2.*fac;
  f3 = -(6*r0.^4-6*r0.^2+1).*fac;
end
if (n >= 4)
  fac = tmpPw.^2.*fac;
  f4 = (24*r0.^6-36*r0.^4+12*r0.^2+1).*fac;
end  
  
if nprime(1)=='0'
  phi = tmp;

elseif nprime(1)=='1'
  % phi = 2ri*f1
  phi = 2.*sq(r(:,:,dim+1)).*f1;

elseif nprime(1)=='2'
  % phi = 4ri^2*f2 + 2*f1
  phi = 4*sq(r(:,:,dim+1)).^2.*f2 + 2*f1;

elseif nprime(1)=='3'
  % phi = 8ri^3*f3+12*ri*f2
  phi = (8*sq(r(:,:,dim+1)).^2.*f3 + 12*f2).*sq(r(:,:,dim+1));
  
elseif nprime(1)=='4'
  % phi = 16*ri^4*f4 + 48ri^2*f3 + 12*f2
  phi = (16*sq(r(:,:,dim+1)).^2.*f4 + 48*f3).*sq(r(:,:,dim+1)).^2 + 12*f2;
  
elseif nprime(1)=='L' && length(nprime)==1
  % phi = 4r^2*f2+2*d*f1
  phi = 4*r0.^2.*f2 + 2*nd*f2;

elseif nprime(1:2)=='L2'
  % phi = 16*r^4*f4 + 16*(d+2)*r^2*f3 + 4d(d+2)*f2
  phi = 16*r0.^4.*f4 + 16*(nd+2)*r0.^2.*f3 + ...
	4*nd*(nd+2)*f2;
  
elseif nprime(1:2)=='m2'
  phi = 4*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).*f2;
  
elseif nprime(1:2)=='m3'
  ndim = length(unique(dim));
  if ndim==2 % dx_i^2 dx_j
      dim = sort(dim);
      if dim(1)==dim(2) % Let d1 be the dimension with a double derivative
          d1 = dim(1); d2 = dim(3);
      else
          d1 = dim(3); d2 = dim(1);
      end
      % 8ri^2rj*f3 + 4rj*f2
      phi = (8*sq(r(:,:,d1+1)).^2.*f3 + 4*f2).*sq(r(:,:,d2+1));
      
  elseif ndim==3 % dx_i dx_j dx_k
      % 8rirjrk*f3
      phi = 8*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).* ...
            sq(r(:,:,dim(3)+1)).*f3;
  else
      error('Error in input argument dim to function bmp')
  end    

elseif nprime(1:2)=='m4'
  dim = sort(dim);
  ndim = length(unique(dim));
  if ndim==2
    if dim(1)==dim(2) & dim(3)==dim(4) % Two derivatives on each coordinate
      % 16ri^2rj^2*f4 + 8*(ri^2+rj^2)*f3 + 4*f2
      phi = 16*sq(r(:,:,dim(1)+1)).^2.*sq(r(:,:,dim(3)+1)).^2.*f4 + ...
	    8*(sq(r(:,:,dim(1)+1)).^2 + sq(r(:,:,dim(3)+1)).^2).*f3 + 4*f2;
          
      else % The 3-1 case
          if dim(1)==dim(2)
              d1 = dim(1); d2 = dim(4);
          else
              d1 = dim(2); d2 = dim(1);
          end
	  % phi = 16ri^3rj*f4+24rirj*f3
          phi = (16*sq(r(:,:,d1+1)).^2.*f4 + 24*f3).* ...
		sq(r(:,:,d1+1)).*sq(r(:,:,d2+1));
          
      end    
  elseif ndim==3 % One coordinate has two derivatives
      if dim(1)==dim(2)
          d1 = dim(1); d2 = dim(3); d3 = dim(4);
      elseif dim(2)==dim(3)
          d1 = dim(2); d2 = dim(1); d3 = dim(4);
      else     
          d1 = dim(3); d2 = dim(1); d3 = dim(2);
      end
      % 16ri^2rjrk*f4+8*rjrk*f3
      phi = (16*sq(r(:,:,d1+1)).^2.*f4 + 8*f3).* ...
	    sq(r(:,:,d2+1)).*sq(r(:,:,d3+1));
      
  elseif ndim==4 % All have one derivative
      % phi = 16rirjrkrl*f4
      phi = 16*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).* ...
	       sq(r(:,:,dim(3)+1)).*sq(r(:,:,dim(4)+1)).*f4;
            
  else
      error('Error in input argument dim to function bump')
  end    
  
  
else
  error('Error in input argument nprime to function bump')
end 
%
% Apply the derivative scaling from rho
%
phi = (1/rho)^n*phi;

function r=sq(r)
 r=squeeze(r);
