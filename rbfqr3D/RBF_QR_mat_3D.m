function [A,T] = RBF_QR_mat_3D(Psi,op,var)
%--- Call patterns
%--- [A,T]=RBF_QR_mat_3D(Psi,op,xe) First time for these points (xe)
%--- [A,T]=RBF_QR_mat_3D(Psi,op,T)  Subsequent times with same xe
%--- [A]=RBF_QR_mat_3D(Psi,op,xe) Just one call with nice interface

if strcmp(op,'0')
    deg = 0;
else
    deg = 1;
end

if isstruct(var)
    T = var;
    xe = T.xe;
else
    T = [];
    xe = var;
end

n = size(xe,1);
N = length(Psi.xk);
j = Psi.j; m = Psi.m; nu = Psi.nu; p = Psi.p; ep = Psi.ep;

T = RBF_QR_precomp_3D(deg,xe,ep,j,m,p,T);
M = length(j);

% jmax = max(j);
% jvec = 0:jmax;
% Tk = cos(acos(xe(:,1))*jvec);
% Yk = Y(jmax,xe(:,2),xe(:,3));
% Pk = ones(n,jmax+1);
% for k=1:jmax
%     Pk(:,k+1) = xe(:,1).*Pk(:,k);
% end
%--- Evaluate the basis functions and compute the interpolation matrix. 
switch(op)
    case '0'
        V = exp(-ep^2*xe(:,1).^2)*ones(1,M).*T.Pk(:,2*m+1).*T.Tk(:,j-2*m+1);
        for k=1:M
            mu = 2*m(k)+p(k);
            if (nu(k)>=0)
                V(:,k) = V(:,k).*T.Yk(:,mu-nu(k)+1,mu+1);
            else
                V(:,k) = V(:,k).*T.Yk(:,mu+1,mu+nu(k)+1);
            end
        end
    case {'x','y','z'}
        V = exp(-ep^2*xe(:,1).^2)*ones(1,M);
        sc3 = zeros(n,2);
        switch(op)
            case 'x'
                sc2 = [sin(xe(:,2)) cos(xe(:,2))];
                sc3(:,1) = cos(xe(:,3));
                sc3(:,2) = sin(xe(:,3));
            case 'y'
                sc2 = [sin(xe(:,2)) cos(xe(:,2))];
                sc3(:,1) = sin(xe(:,3));
                sc3(:,2) = -cos(xe(:,3));
            case 'z'
                sc2 = [cos(xe(:,2)) -sin(xe(:,2))];
                sc3(:,1) = ones(n,1);
        end
        for k=1:M
            mk = m(k);
            mu = 2*mk+p(k);
            if nu(k)>=0
                ii = mu-nu(k)+1; jj = mu+1;
            else
                ii = mu+1; jj = mu+nu(k)+1;
            end
            if mk>0
                V(:,k) = V(:,k).*( ...
                  T.f_prime(:,k).*sc3(:,1).*sc2(:,1).*T.Yk(:,ii,jj) - ...
                  T.f_over_r(:,k).*sc3(:,2).*T.Yfi_over_s(:,ii,jj) + ...
                  T.f_over_r(:,k).*sc3(:,1).*sc2(:,2).*T.Yth(:,ii,jj));
            else
                if p(k)==0
                  V(:,k) = V(:,k).*( ...
                    T.f_prime(:,k).*sc3(:,1).*sc2(:,1).*T.Yk(:,ii,jj));
                else
                  V(:,k) = V(:,k).*( ...
                    T.f_prime(:,k).*sc3(:,1).*sc2(:,1).*T.Yk(:,ii,jj) -...
                    T.f_over_r(:,k).*sc3(:,2).*T.Yfi_over_s(:,ii,jj) + ...
                    T.f_over_r(:,k).*sc3(:,1).*sc2(:,2).*T.Yth(:,ii,jj));
                end
            end
        end
    case {'xx','yy','zz'}
        V = exp(-ep^2*xe(:,1).^2)*ones(1,M);
        % r, th, fi
        cosfi = cos(xe(:,3)); sinfi = sin(xe(:,3));
        costh = cos(xe(:,2)); sinth = sin(xe(:,2));
        switch(op)
          case 'xx'
            % Multiplying no derivative, one deriv, and two derivs
            sc0 = [cosfi.^2.*sinth.^2 sinfi.^2+cosfi.^2.*costh.^2];
            sc1 = [sinfi.*cosfi cosfi.^2.*sinth.*costh];
            sc2 = [cosfi.^2.*costh.^2 sinfi.^2 sc1(:,1)];
          case 'yy'
            sc0 = [sinfi.^2.*sinth.^2 cosfi.^2+sinfi.^2.*costh.^2];
            sc1 = [-sinfi.*cosfi sinfi.^2.*sinth.*costh];
            sc2 = [sinfi.^2.*costh.^2 cosfi.^2 sc1(:,1)];
          case 'zz'
            sc0 = [costh.^2 sinth.^2];
            sc1 = zeros(n,2);
            sc1(:,2) = -sinth.*costh;
            sc2 = zeros(n,3);
            sc2(:,1) = sinth.^2;
        end
 
        for k=1:M
          mk = m(k);
          mu = 2*mk+p(k);
          if nu(k)>=0
            ii = mu-nu(k)+1; jj = mu+1;
          else
            ii = mu+1; jj = mu+nu(k)+1;
          end
          if mk>0            
            V(:,k) =  V(:,k).*(...
              (sc0(:,1).*T.f_bis(:,k) + sc0(:,2).*T.f_prime_over_r(:,k)).* ...
                T.Yk(:,ii,jj) ...
              -2*sc1(:,1).*T.f_prime_over_r(:,k).*T.Yfi(:,ii,jj) ...
              +2*sc1(:,2).*(T.f_prime_over_r(:,k)-T.f_over_r2(:,k)).* ...
                T.Yth(:,ii,jj) ...
              +sc2(:,1).*T.f_over_r2(:,k).*T.Ythth(:,ii,jj) ...
              +sc2(:,2).*T.f_over_r2(:,k).*T.Yfifi(:,ii,jj) ...
              +2*sc2(:,3).*T.f_over_r2(:,k).*T.Ythfi(:,ii,jj));
          else
            if p(k)==0
              V(:,k) = V(:,k).*( ...
                (sc0(:,1).*T.f_bis(:,k)+sc0(:,2).*T.f_prime_over_r(:,k)).* ...
                T.Yk(:,ii,jj));
            else
              V(:,k) =  V(:,k).*(...
              (sc0(:,1).*T.f_bis(:,k) ...
               + sc0(:,2).*T.f_prime_over_r_minus_f_over_r2(:,k)).* ...
                T.Yk(:,ii,jj) ...
              +T.f_prime_over_r_minus_f_over_r2(:,k).* ...
                 (-2*sc1(:,1).*T.Yfi(:,ii,jj) + 2*sc1(:,2).*T.Yth(:,ii,jj)));
            end
          end  
        end
  case {'xy','yx','zx','xz','zy','yz'}
        V = exp(-ep^2*xe(:,1).^2)*ones(1,M);
        % r, th, fi
        cosfi = cos(xe(:,3)); sinfi = sin(xe(:,3));
        costh = cos(xe(:,2)); sinth = sin(xe(:,2));
        % Cofficients for f''Y, (f'/r)Y
        % (f'/r)Yfi, (f'/r)Yfi/sinth, (f'/r-f/r^2)Yth
        % (f/r^2)Ythh, (f/r^2)Yfifi, (f/r^2)Ythfi, real Ythfi
        switch(op)
          case {'xy','yx'}
            sc0 = [(sinfi.*cosfi.*sinth.^2) (-sinfi.*cosfi.*sinth.^2)];
            sc1 = [(cosfi.^2-sinfi.^2) zeros(n,1) ...
                   2*(sinfi.*cosfi.*sinth.*costh)];
            sc2 = [(sinfi.*cosfi.*costh.^2) (-sinfi.*cosfi) ...
                   -(cosfi.^2-sinfi.^2) zeros(n,1)];
          case {'xz', 'zx'}
            sc0 = [(cosfi.*sinth.*costh) -(cosfi.*sinth.*costh)];
            sc1 = [zeros(n,1) (-sinfi.*costh) cosfi.*(costh.^2-sinth.^2) ];
            sc2 = [-(cosfi.*sinth.*costh) zeros(n,1) zeros(n,1) sinfi];
          case {'yz', 'zy'}
            sc0 = [(sinfi.*sinth.*costh) -(sinfi.*costh.*sinth)];
            sc1 = [zeros(n,1) (cosfi.*costh) sinfi.*(costh.^2-sinth.^2)];
            sc2 = [-(sinfi.*sinth.*costh) zeros(n,1) zeros(n,1) -cosfi];
        end
 
        for k=1:M
          mk = m(k);
          mu = 2*mk+p(k);
          if nu(k)>=0
            ii = mu-nu(k)+1; jj = mu+1;
          else
            ii = mu+1; jj = mu+nu(k)+1;
          end
          if mk>0            
            V(:,k) =  V(:,k).*(...
              (sc0(:,1).*T.f_bis(:,k) + sc0(:,2).*T.f_prime_over_r(:,k)).* ...
                T.Yk(:,ii,jj) ...
              + T.f_prime_over_r(:,k).* ...
              (sc1(:,1).*T.Yfi(:,ii,jj) + sc1(:,2).*T.Yfi_over_s(:,ii,jj)) ...
              +sc1(:,3).*(T.f_prime_over_r(:,k)-T.f_over_r2(:,k)).* ...
                T.Yth(:,ii,jj) ...
              +T.f_over_r2(:,k).* ...
            (sc2(:,1).*T.Ythth(:,ii,jj) + sc2(:,2).*T.Yfifi(:,ii,jj) ...
             +sc2(:,3).*T.Ythfi(:,ii,jj)) ...
              +T.f_over_r2(:,k).*sc2(:,4).*(-nu(k)).*T.Yth(:,jj,ii));
          else
            if p(k)==0
              V(:,k) = V(:,k).*( ...
                (sc0(:,1).*T.f_bis(:,k) + sc0(:,2).*T.f_prime_over_r(:,k)).* ...
                T.Yk(:,ii,jj));
            else
              V(:,k) =  V(:,k).*(...
                (sc0(:,1).*T.f_bis(:,k) + sc0(:,2).* ...
                 T.f_prime_over_r_minus_f_over_r2(:,k)).*T.Yk(:,ii,jj) ...
                +T.f_prime_over_r_minus_f_over_r2(:,k).* ...
                (sc1(:,1).*T.Yfi(:,ii,jj) + sc1(:,2).*T.Yfi_over_s(:,ii,jj)) ...
                +sc1(:,3).*T.f_prime_over_r_minus_f_over_r2(:,k).* ...
                T.Yth(:,ii,jj));
            end
          end  
        end

    case 'L'
        V = exp(-ep^2*xe(:,1).^2)*ones(1,M);
        for k=1:M
            mk = m(k);
            mu = 2*mk+p(k);
            if nu(k)>=0
                ii = mu-nu(k)+1; jj = mu+1;
            else
                ii = mu+1; jj = mu+nu(k)+1;
            end
            if mk>0
                V(:,k) = V(:,k).*(T.f_bis(:,k) + ...
                             2*T.f_prime_over_r(:,k) - mu*(mu+1)*T.f_over_r2(:,k)).*...
                             T.Yk(:,ii,jj);                
            else
                if p(k)==0
                            V(:,k) = V(:,k).*(T.f_bis(:,k)+2*T.f_prime_over_r(:,k)).*T.Yk(:,ii,jj);
                else
                            V(:,k) = V(:,k).*(T.f_bis(:,k)+2*T.f_prime_over_r_minus_f_over_r2(:,k)).*...
                                T.Yk(:,ii,jj);
                end
            end
        end
    case {'L1','L2','L3','L4'}
        V = zeros(n,M);
        d = str2num(op(2));
        Lpow_coef = load('Lpow_coef_j40_k4.mat','T'); % Load saved coefficients
        for k=1:M
            jk = j(k);
            mk = m(k);
            mu = 2*mk+p(k);
            if nu(k)>=0
                ii = mu-nu(k)+1; jj = mu+1;
            else
                ii = mu+1; jj = mu+nu(k)+1;
            end
            fe = Lpow_eval(Lpow_coef.T,xe(:,1),ep,jk,mk,d);
            V(:,k) = fe.*T.Yk(:,ii,jj);
        end
end
V = V(:,Psi.columns);
A = V(:,1:N) + V(:,N+1:end)*Psi.Rt';
