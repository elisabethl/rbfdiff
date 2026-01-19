function [phi]=w2(rho,r,nprime,dim)
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
%
% Only derivatives of order at most 2 are allowed for this Wendland C^2 function.
%
switch nprime
  case {'0','1','2'}
    n = str2num(nprime);				
  case {'L','m2'}
    n = 2;
  case {'3','4','L2','m3','m4'}
    error('Derivatives of order higher than 2 are not allowed for w2')
  otherwise
    error('Unknown operator nprime in phs')
end
  
if nargin<=3
  dim=1;
end
%
% For the case of the L operator, we need to know the number of dimensions
%
if (nprime(1)=='L')
  if (size(r,3)==1)
    nd=dim;
  else
    nd=size(r,3)-1;
  end
end
%
% For the mixed derivative, the dimensions must be given even in 2D
%
if (nprime(1)=='m')
  if (length(dim)~=2)
    error('For the mixed derivative, dim=dim(1:2)')
  elseif (dim(1)==dim(2))
    error('For mixed derivatives, dim(1) must be other than dim(2)')
  end  
end

phi=zeros(size(r,1),size(r,2));
%
% Here we scale the given radius with rho, such that we can compute all derivatives in the reference scaling. Later we add the scaling from differentiation.
r = (1/rho)*r;

%
% Extract the parts of r that will be used
%
r0 = sq(r(:,:,1));
sz = size(r0);
r0 = r0(:);
pos0 = find(r0 <= 0);
if nprime(1)~='L' & nprime(1)~='0' % Else only radial powers are needed
    for d=1:length(dim)
        ri = sq(r(:,:,dim(d)+1));
        riMat(:,dim(d)) = ri(:);
        riDiv(:,dim(d)) = ri(:)./r0;
        riDiv(pos0,dim(d)) = 1; % Limit value coming from the positive side
    end
end
clear r
%
% Make sure the result is zero for any value with radius larger than 1.
%
tmp = (1-r0).*(r0 <= 1);

if nprime(1)=='0'
  phi = (tmp.^4).*(4*r0 + 1);

elseif nprime(1)=='1'
  % -20*ri*tmp^3
  phi = -20*riMat(:,dim).*(tmp.^3);

elseif nprime(1)=='2'
  % 4ri^2*f2 + 2*f1
  phi = 60*riDiv(:,dim).*riMat(:,dim).*tmp.^2 - 20*tmp.^3;

elseif nprime(1)=='L' && length(nprime)==1
  % 4r^2*f2+2*d*f1
  phi = 60*r0.*tmp.^2 + -20*nd*tmp.^3;

elseif nprime(1:2)=='m2'
  % 4ri*rj*f2
  phi = 60*riDiv(:,dim(1)).*riDiv(:,dim(2)).*r0.*tmp.^2;

else
  error('Error in input argument nprime to function w2')
end

phi = reshape(phi,sz);

phi = (1/rho)^n*phi;

function r=sq(r)
 r=squeeze(r);
