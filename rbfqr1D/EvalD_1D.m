function D=EvalD_1D(ep,p1,p2,j)
% The purpose of this function is to compute the scaling effect of 
% D_1^{-1} and D_2 applied to the correction matrix in the RBF-QR method.
%--- ep (scalar) : The shape parameter
%--- p1 (vector) : Indices for elements in D_1 
%--- p2 (vector) : Indices for elements in D_2
%--- j  (vector) : Identifiers for the expansion functions T_{j}  
%--- ep^{2*(j2-j1)}*j1!/j2!
   
%--- An assumption here, associated with the RBF-QR method is that j_2>=j_1  
  D = zeros(length(p1),length(p2));
  
%--- Precompute the powers of epsilon that will be needed 
  ep2 = ep*ep;
  epp(1) = 1; %ep^0
  for pp=1:j(p2(end))
    epp(pp+1) = ep2*epp(pp);
  end  
%--- Put the powers into D row by row. The column powers can be reused.
  pow2 = j(p2);
  for k=1:length(p1)
    D(k,:) = epp(pow2-j(p1(k)) +1);
  end
%--- Then we do the factorials using precomputed values and dividing  
  fp = cumprod([1 1:j(p2(end))]);
  denom = 1./fp(j(p2) +1); % The inverse denominator actually
  for k=1:length(p1)
    D(k,:) = D(k,:).*(denom*fp(j(p1(k)) +1));
  end
  