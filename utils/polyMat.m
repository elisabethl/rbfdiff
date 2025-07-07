% -------------------------------------------------------------------------
% polyMat
% Purpose: Given points, a polynomial degree, and a vector with the number
%          of derivatives in each dimension, compute a monomial basis matrix
%          evaluated at these points with the derivatives applied to each
%          monimial, represented by one column.
% P = polyMat(xc,pdeg,diff)
% Input:  xc     Double N x d array of N points in d dimensions 
% 
%         pdeg	 Integer scalar, the maximum polynomial degree
%         diff   Integer 1 x d vector of derivative orders
%
% Output: P	 Double N x dim(P,d) matrix of the potentially differentiated
%                monomial matrix avaluated at xc.
% Example:       For 10 data points in 3D, compute the first derivative in the
%                second coordinate of all monomials of degree up to 5.
%                xc = rand(10,3); 
%                pdeg = 5; 
%                diff = [0 1 0];
%                P = polyMat(xc,pdeg,diff);
% Copyright (c) 2024 Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%		       Andreas Michael <andreas.michael@it.uu.se >
% -------------------------------------------------------------------------
function P = polyMat(xc,pdeg,diff)
  if (pdeg>=0) % Note that we allow calls with, e.g.,  pdeg = -1
    nd = size(xc,2);
    nc = size(xc,1);
    %
    % First, we get the powers for the columns of the matrix
    %
    pow = polyPow(pdeg,nd);
    %
    % Then we construct base matrices for each dimension
    %
    pmat = zeros(nc,nd,pdeg+1);
    pmat(:,:,1) = 1;
    for k=2:pdeg+1
      pmat(:,:,k) = pmat(:,:,k-1).*xc;
    end
    %
    % Then build P
    %
    P = ones(nc,size(pow,1));
    for d=1:nd
      if (nargin==3)
        if(diff(d)>0)
          %
          % Recompute the powers according to derivatives, and compute the coeff
          %
          coeff = pow(:,d);
          for m=2:diff(d)
            % If we reach zero, the coefficient becomes zero after that
            coeff = coeff.*(pow(:,d)-m+1);
          end  
          pow(:,d) = pow(:,d) - diff(d);
          pow(:,d) = max(0,pow(:,d)); % just to avoid computing negative powers
          % Multiply with coefficient due to differentiation in dim=d 
          P = P.*repmat(coeff(:)',nc,1);
        end
      end
      % Multiply each column with the appropriate power of x_d
      P = P.*reshape(pmat(:,d,pow(:,d)+1),size(P));
    end  
  else
    P = zeros(size(xc,1),0);
  end
  
  
