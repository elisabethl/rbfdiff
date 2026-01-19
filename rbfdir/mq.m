function [phi]=mq(epsil,r,nprime,dim)
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
% For the mixed derivatives, the dimensions must be given for all d
%
if (nprime(1)=='m')
  p = str2num(nprime(2));  
  if (length(dim)~=p)
    error('For the mixed derivative mp, dim=dim(1:p)')
  elseif (all(dim==dim(1)))
    error('For mixed derivatives, dim(:) cannot have only one value')
  end  
end
%
% epsil can be either just one value or a vector of N values
%
esz = size(epsil);
if (prod(esz)~=1)
  if (min(esz)==1 & max(esz)==size(r,2))
    %
    % Make epsil into a matrix with constant columns
    %
    epsil = ones(size(r,1),1)*epsil(:).';
  else
    error('The size of epsil does not match the columns in r')
  end
end

phi=zeros(size(r,1),size(r,2));

tmp = 1+(epsil.*sq(r(:,:,1))).^2;

if nprime(1)=='0'
  phi = sqrt(tmp);

elseif nprime(1)=='1'
  phi = epsil.^2.*sq(r(:,:,dim+1))./sqrt(tmp);

elseif nprime(1)=='2'
  phi = epsil.^2./sqrt(tmp) - epsil.^4.*sq(r(:,:,dim+1)).^2.*tmp.^-1.5;

elseif nprime(1)=='3'
    phi = -3*epsil.^4.*sq(r(:,:,dim+1))   .*tmp.^-1.5 + ...
           3*epsil.^6.*sq(r(:,:,dim+1)).^3.*tmp.^-2.5;

elseif nprime(1)=='4'
  phi = -3*epsil.^4.*tmp.^-1.5 +18*epsil.^6.*sq(r(:,:,dim+1)).^2.*tmp.^-2.5 - ...
                              15*epsil.^8.*sq(r(:,:,dim+1)).^4.*tmp.^-3.5;

elseif nprime(1)=='L' & length(nprime)==1
  phi = nd*epsil.^2.*tmp.^-0.5 - epsil.^4.*sq(r(:,:,1)).^2.*tmp.^-1.5;

elseif nprime(1:2)=='L2'
  phi = -nd*(nd+2)*epsil.^4.*tmp.^-1.5 + ...
         6*(nd+2)*epsil.^6.*sq(r(:,:,1)).^2.*tmp.^-2.5 - ...
         15*epsil.^8.*sq(r(:,:,1)).^4.*tmp.^-3.5; 

elseif nprime(1:2)=='m2'
  phi = -epsil.^4.*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).*tmp.^-1.5; 

elseif nprime(1:2)=='m3'
  ndim = length(unique(dim));
  if ndim==2 % dx_i^2 dx_j
      dim = sort(dim);
      if dim(1)==dim(2) % Let d1 be the dimension with a double derivative
          d1 = dim(1); d2 = dim(3);
      else
          d1 = dim(3); d2 = dim(1);
      end    
      phi = 3*epsil.^6.*sq(r(:,:,d1+1)).^2.*sq(r(:,:,d2+1)).*tmp.^-2.5 - ...
              epsil.^4.*sq(r(:,:,d2+1)).*tmp.^-1.5;
      
  elseif ndim==3 % dx_i dx_j dx_k
      phi = 3*epsil.^6.*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).* ...
            sq(r(:,:,dim(3)+1)).*tmp.^-2.5;
  else
      error('Error in input argument dim to function mq')
  end    

elseif nprime(1:2)=='m4'
  dim = sort(dim);
  ndim = length(unique(dim));
  if ndim==2
      if dim(1)==dim(2) & dim(3)==dim(4) % Two derivatives on each coordinate
          phi = -15*epsil.^8.*sq(r(:,:,dim(1)+1)).^2.* ...
                              sq(r(:,:,dim(3)+1)).^2.*tmp.^-3.5 + ...
                  3*epsil.^6.*(sq(r(:,:,dim(1)+1)).^2 + ...
                               sq(r(:,:,dim(3)+1)).^2).*tmp.^-2.5 -  ...
                    epsil.^4.*tmp.^-1.5;
          
      else % The 3-1 case
          if dim(1)==dim(2)
              d1 = dim(1); d2 = dim(4);
          else
              d1 = dim(2); d2 = dim(1);
          end
          phi = -15*epsil.^8.*sq(r(:,:,d1+1)).^3.*sq(r(:,:,d2+1)).*tmp.^-3.5 + ...
                 9*epsil.^6.*sq(r(:,:,d1+1)).*sq(r(:,:,d2+1)).*tmp.^-2.5;
          
      end    
  elseif ndim==3 % One coordinate has two derivatives
      if dim(1)==dim(2)
          d1 = dim(1); d2 = dim(3); d3 = dim(4);
      elseif dim(2)==dim(3)
          d1 = dim(2); d2 = dim(1); d3 = dim(4);
      else     
          d1 = dim(3); d2 = dim(1); d3 = dim(2);
      end
      phi = -15*epsil.^8.*sq(r(:,:,d1+1)).^2.*sq(r(:,:,d2+1)).* ...
                          sq(r(:,:,d3+1)).*tmp.^-3.5 +  ...
              3*epsil.^6.*sq(r(:,:,d2+1)).*sq(r(:,:,d3+1)).*tmp.^-2.5;
      
  elseif ndim==4 % All have one derivative
      phi = -15*epsil.^8.*sq(r(:,:,dim(1)+1)).*sq(r(:,:,dim(2)+1)).* ...
                          sq(r(:,:,dim(3)+1)).*sq(r(:,:,dim(4)+1)).*tmp.^-3.5;
  else
      error('Error in input argument dim to function mq')
  end    
  
else
  error('Error in input argument nprime to function mq')
end 

function r=sq(r)
 r=squeeze(r);
