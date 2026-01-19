function [phi]=phs(m,r,nprime,dim)
% This is a polyharmonic spline (PHS) phi(r) = r^m.
% NPRIME is a string defining which operator to use on the basis function
%
% DIM is the dimension for the partial derivative if nprime is '0','1'...,'4'
% DIM(1:2) are the dimensions for the mixed second derivative if nprime
% is 'm2'
% and DIM is the number of space dimensions if nprime is 'L' or 'L2'
%
% Assume a one-dimensional problem if no dimension is given.
% Copyright 2016, Slobodan Milovanovic based on Elisabeth Larsson.
% Modified with more derivatives 2026, Elisabeth Larsson

%
% Check that the order is an odd number
%
if mod(m,2)~=1
    m
    error('Only odd powers for r^m are accepted in the phs subroutine')
end
%
% Check that the derivatives exist. We will allow also the discontinuous
% derivative that can be taken one time after the continuous ones.
%
switch nprime
  case {'0','1','2','3','4'}
    n = str2num(nprime);
  case 'L'
    n = 2;
  case 'L2'
    n = 4;
  case {'m2','m3','m4'}
    n = str2num(nprime(2));
  otherwise
    error('Unknown operator nprime in phs')
end
if n>m
    error('No more than m derivatives can be taken of the phs basis functions')
end

if nargin <= 3
    dim = 1;
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
    p = str2num(nprime(2));
    if (length(dim)~=p)
        error('For the mixed derivative mp, dim=dim(1:p)')
    elseif (all(dim==dim(1)))
        error('For mixed derivatives, dim(:) cannot have only one value')
    end
end
%
% Use ri/r to avoid singularity. We extract the different powers.
%
r0 = sq(r(:,:,1));
sz = size(r0);
r0 = r0(:);
pos0 = find(r0<=0); % All values should be positive, but in case, we take also <0.

if nprime(1)~='L' & nprime(1)~='0' % Else only radial powers are needed
    dim = sort(dim);
    udim = unique(dim);
    for d=1:length(udim)
        ri = sq(r(:,:,udim(d)+1));
        riDiv(:,udim(d)) = ri(:)./r0;
        riDiv(pos0,udim(d)) = 1; % Limit value coming from the positive side
    end
end
clear r

if nprime(1)=='0'
    phi = r0.^m;
    
elseif nprime(1)=='1'
    % phi = ri/r*m*r^(m-1)
    phi = m*r0.^(m-1).*riDiv(:,dim);
    
elseif nprime(1)=='2'
    % phi = (ri/r)^2m*(m-2)r^(m-2) + m*r^(m-2)
    phi = m*(m-2)*r0.^(m-2).*riDiv(:,dim).^2 + m*r0.^(m-2);

elseif nprime(1)=='3'
    % phi = m*(m-2)*(m-4)r^(m-3)(ri/r)^3 + 3*m*(m-2)r^(m-3)(ri/r)
    phi = m*(m-2)*((m-4)*riDiv(:,dim).^3 + 3*riDiv(:,dim)).*r0.^(m-3);

elseif nprime(1)=='4'
    % phi = ri^4*m*(m-2)*(m-4)*(m-6)*r^(m-8) + 6*ri^2*m*(m-2)*(m-4)r^(m-6)
    % + 3*m*(m-2)*r^(m-4) 
    phi = m*(m-2)*((m-4)*(m-6)*riDiv(:,dim).^4 + 6*(m-4)*riDiv(:,dim).^2 + 3)*r0.^(m-4);

elseif nprime(1)=='L' & length(nprime)==1
  % phi = m*(m-2)r^{m-2} + d*m*r^(m-2)
  phi = m*(m-2+nd).*r0.^(m-2);

elseif nprime(1:2)=='L2'
  % phi = m*(m-2)*(m-4)*(m-6)r^(m-4) + 2(d+2)m*(m-2)*(m-4)*r^(m-4) +
  % d(d+2)m*(m-2)*r^(m-4)
  phi = m*(m-2)*((m-4)*(m-6)+2*(nd+2)*(m-4)+nd*(nd+2))*r0.^(m-4); 
    
elseif nprime(1:2)=='m2'
  % phi = r1*rj*m*(m-2)r^{m-4}
  phi = m*(m-d)*riDiv(:,dim(1)).*riDiv(:,dim(2)).*r0.^(m-2);

elseif nprime(1:2)=='m3'
  ndim = length(unique(dim));
  if ndim==2 % dx_i^2 dx_j
      dim = sort(dim);
      if dim(1)==dim(2) % Let d1 be the dimension with a double derivative
          d1 = dim(1); d2 = dim(3);
      else
          d1 = dim(3); d2 = dim(1);
      end
      % phi = m*(m-2)*(m-4)ri^2rjr^{m-6}+ m*(m-2)rj*r^(m-4)
      phi = m*(m-2)*riDiv(:,d2).*((m-4)*riDiv(:,d1).^2 + 1).*r0.^(m-3);
      
  elseif ndim==3 % dx_i dx_j dx_k
      % phi = rirjrk m*(m-2)*(m-4)*r^(m-6)
    phi = m*(m-2)*(m-4)*riDiv(:,dim(1)).*riDiv(:,dim(2)).*riDiv(:,dim(3)).*r0.^(m-3);
  else
      error('Error in input argument dim to function mq')
  end    

elseif nprime(1:2)=='m4'
  dim = sort(dim);
  ndim = length(unique(dim));
  if ndim==2
    if dim(1)==dim(2) & dim(3)==dim(4) % Two derivatives on each coordinate
        % phi = 16ri^2rj^2*f4 + 8(ri^2+rj^2)*f3 + 4*f2
      phi = m*(m-2)*((m-4)*(m-6)*riDiv(:,dim(1)).^2.*ridiv(:,dim(3)).^2 + ...
	             (m-4)*(riDiv(:,dim(1)).^2 + riDiv(:,dim(3)).^2) + 1).*r0.^(m-4);
                        
      else % The 3-1 case
          if dim(1)==dim(2)
              d1 = dim(1); d2 = dim(4);
          else
              d1 = dim(2); d2 = dim(1);
          end
	  % phi = 16ri^3rj*f4 + 24ri*rj*f3
          phi = m*(m-2)*(m-4)*riDiv(:,d1).*riDiv(:,d2).*( ...
		(m-6)*riDiv(:,d1).^2 + 3).*r0.^(m-4); 
      end    
  elseif ndim==3 % One coordinate has two derivatives
      if dim(1)==dim(2)
          d1 = dim(1); d2 = dim(3); d3 = dim(4);
      elseif dim(2)==dim(3)
          d1 = dim(2); d2 = dim(1); d3 = dim(4);
      else     
          d1 = dim(3); d2 = dim(1); d3 = dim(2);
      end
      % phi = 16ri^2rjrk*f4 + 8rj*rk*f3
      phi = m*(m-2)*(m-4)*riDiv(:,d2).*riDiv(:,d3).*((m-6)*riDiv(:,d1).^2 + 1).*r0.^(m-4);
            
  elseif ndim==4 % All have one derivative
      % phi = 16rirjrkrl*f4
      phi = m*(m-2)*(m-4)*(m-6)*riDiv(:,dim(1)).*riDiv(:,dim(2)).* ...
	                        riDiv(:,dim(3)).*riDiv(:,dim(4)).*r0.^(m-4);
            
  else
      error('Error in input argument dim to function mq')
  end        
    
else
    error('Error in input argument nprime to function phs')
end

phi = reshape(phi,sz);

function r=sq(r)
r=squeeze(r);
