function [phi]=wendland_c2(epsil,r,nprime,dim)
%
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% DIM(1:2) are the dimensions for the mixed second derivative if nprime
% is 'm2'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
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

tmp = (1-sq(r(:,:,1)));

if nprime(1)=='0'
  phi = (tmp.^4).*(4.*sq(r(:,:,1)) + 1).*(sq(r(:,:,1))<=1);

elseif nprime(1)=='1'
  phi = -20.*sq(r(:,:,dim+1)).*(tmp.^3).*(sq(r(:,:,1))<=1);

elseif nprime(1)=='2'
  % Making sure x_i^2/r = 0 at r = 0.
  phi = (60.*(sq(r(:,:,dim+1).^2)./(sq(r(:,:,1)) + (sq(r(:,:,1)) == 0))).*tmp.^2 - 20.*tmp.^3).*(sq(r(:,:,1))<=1);

elseif nprime(1)=='L' && length(nprime)==1
  phi = (60.*sq(r(:,:,1)).*tmp.^2 - nd.*20.*tmp.^3).*(sq(r(:,:,1))<=1);

elseif nprime(1:2)=='m2'
  % Making sure x_ix_j/r = 0 at r = 0.
  phi = (60.*(sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1))./(sq(r(:,:,1)) + (sq(r(:,:,1)) == 0))).*tmp.^2).*(sq(r(:,:,1))<=1);
else
  error('Error in input argument nprime to function bump')
end 

function r=sq(r)
 r=squeeze(r);
