% -------------------------------------------------------------------------
% polyPow
% Purpose: Returns all the monomial powers involved in a polynomial of degree
%          pdeg in d dimensions.
% pow = polyPow(pdeg,d)
% Input:  pdeg	 Integer scalar, the maximum polynomial degree
%         d      Integer scalar, the number of dimensions
%
% Output: pow	 integer N x d matrix of the powers of the N monomials
%
% Example:       pdeg = 4;
%                d = 3;
%                pow = polyPow(pdeg,d);
% Copyright (c) 2024 Elisabeth Larsson <elisabeth.larsson@it.uu.se>
%		       Andreas Michael <andreas.michael@it.uu.se >
% -------------------------------------------------------------------------
function pow = polyPow(p,d)
  for k=1:p+1
    pp{1,k} = k-1;
    n(1,k) = 1;
    for j=2:d
      pp{j,k} = [pp{j-1,k} zeros(n(j-1,k),1)]; 
      for m=k-1:-1:1
        pp{j,k} = [pp{j,k}; pp{j-1,m} (k-m)*ones(n(j-1,m),1)];
      end  
      n(j,k) = size(pp{j,k},1);
    end  
  end
  pow = pp{d,1};
  for k=2:p+1
    pow = [pow;pp{d,k}];  
  end
  
