function [A,T]=RBF_QR_mat_1D(Psi,op,var)
% Compute an RBF-QR matrix for the operator op, that is a_ij=L Psi_j(xe(i))
%--- op (string) : Can be '0', 'x', or 'xx' (which derivative)
%--- xe(1:Ne) : The evaluation points 

%--- Find out which call sign we got
if (isstruct(var))
  T=var;
  xe = T.xe;
else
  T=[];
  xe=var;
end

%--- Psi contains all the information about the RBF-QR basis    
  Ne = size(xe,1);
  ep = Psi.ep;
  Rt = Psi.Rt;
  jmax = Psi.jmax;
  M = Psi.M;
  N = size(Rt,1);

%--- Transform xe to 'polar' coordinates
  
  re = abs(xe);
  se = ones(size(xe));
  pos = find(xe<0);
  se(pos) = -1;

%--- Evaluate the Chebyshev polynomials and their derivatives stably
  if (isempty(T))
      T.xe = xe;
      [T.Te,T.dTe,T.d2Te] = ClenshawCheb(re,jmax);
  
%--- These are only needed for the derivatives, but for clarity, compute here
      T.evec = exp(-ep^2*re.^2);
      T.evecx = 2*ep^2*re.*T.evec;
      T.evecx2 = 4*ep^4*re.^2.*T.evec;
  end

  if (strcmp(op,'0'))
      Ve = (T.evec*ones(1,M)).*T.Te(:,1:M);
      % Multiply by the sign for all odd j
      for k=2:2:jmax+1
          Ve(:,k) = se.*Ve(:,k);
      end
    
  elseif (strcmp(op,'x'))
     
      Ve = (T.evec*ones(1,M)).*T.dTe(:,1:M)-(T.evecx*ones(1,M)).*T.Te(:,1:M);

      for k=1:2:jmax+1
          % Inner derivative dr/dx=s^{-1} combined with s^p
          Ve(:,k) = se.*Ve(:,k);
      end
	
  elseif (strcmp(op,'xx'))
      
      Ve = (T.evec*ones(1,M)).*T.d2Te(:,1:M)- ...
           (2*T.evecx*ones(1,M)).*T.dTe(:,1:M)+ ...
           ((T.evecx2-2*ep^2*T.evec)*ones(1,M)).*T.Te(:,1:M);
      
      for k=2:2:jmax+1
          Ve(:,k) = se.*Ve(:,k);
      end	
  else
        error(['Operator ' op ' not implemented'])
  end    

%--- Evaluate the RBF-QR matrix for the given operator 
  A = Ve(:,1:N) + Ve(:,N+1:M)*Rt.';


