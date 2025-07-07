function T = RBF_QR_precomp_3D(deg,xe,ep,j,m,p,T)
% Precompute some useful things needed for derivatives. Store in T.
jmax = max(j);
jvec = 0:jmax;
n = size(xe,1);

if nargin<7 || isempty(T)
    % Always compute degree zero variables, unless already present.
    T.xe = xe;
    T.Tk = cos(acos(xe(:,1))*jvec);
    
    [T.Yk,T.Yth,T.Yfi_over_s,T.Yfi,T.Ythth,T.Ythfi,T.Yfifi] =  ...
        Y(jmax,xe(:,2),xe(:,3));
    
    T.Pk = ones(n,jmax+3);
    for k=1:jmax+2
        T.Pk(:,k+1) = xe(:,1).*T.Pk(:,k);
    end
    T.deg = 0;
end

if T.deg<deg
    % Currently only degree zero or degree > 0
    rn1 = find(xe(:,1)+eps<1);
    T.deg = deg;
    T.Tk_p = ones(n,1)*jvec.^2;
    T.Tk_p(rn1,:) = repmat(jvec,length(rn1),1).*sin(acos(xe(rn1,1))*jvec)...
        ./repmat(sin(acos(xe(rn1,1))),1,jmax+1);
    T.Tk_b = ones(n,1)*(jvec.^2.*(jvec.^2-1)./3);
    T.Tk_b(rn1,:) = (repmat(jvec.^2,length(rn1),1).*T.Tk(rn1,:)-...
        repmat(xe(rn1,1),1,jmax+1).*T.Tk_p(rn1,:))./...
        repmat(xe(rn1,1).^2-1,1,jmax+1);
    
    M = length(m);
    T.f_over_r = zeros(n,M); T.f_over_r2 = zeros(n,M); T.f_prime = zeros(n,M);
    T.f_prime_over_r = zeros(n,M); T.f_bis = zeros(n,M);
    T.f_prime_over_r_minus_f_over_r2 = zeros(n,M);
    for k=1:M
        mk = m(k);
        jk = j(k);
        if mk>0
            T.f_over_r(:,k) = T.Pk(:,2*mk-1+1).*T.Tk(:,jk-2*mk+1);
            T.f_over_r2(:,k) = T.Pk(:,2*mk-2+1).*T.Tk(:,jk-2*mk+1);
            T.f_prime(:,k) = 2*(mk*T.Pk(:,2*mk-1+1)-ep^2*T.Pk(:,2*mk+1+1)).*T.Tk(:,jk-2*mk+1)+...
                T.Pk(:,2*mk+1).*T.Tk_p(:,jk-2*mk+1);
            T.f_prime_over_r(:,k) = 2*(mk*T.Pk(:,2*mk-2+1)-ep^2*T.Pk(:,2*mk+1)).*T.Tk(:,jk-2*mk+1)+...
                T.Pk(:,2*mk-1+1).*T.Tk_p(:,jk-2*mk+1);
            T.f_bis(:,k) = 2*(2*ep^4*T.Pk(:,2*mk+2+1)-(4*mk+1)*ep^2*T.Pk(:,2*mk+1)+...
                mk*(2*mk-1)*T.Pk(:,2*mk-2+1)).*T.Tk(:,jk-2*mk+1) + ...
                4*(mk.*T.Pk(:,2*mk-1+1)-ep^2*T.Pk(:,2*mk+1+1)).*T.Tk_p(:,jk-2*mk+1) + ...
                T.Pk(:,2*mk+1).*T.Tk_b(:,jk-2*mk+1);
        else
            T.f_prime(:,k) = T.Tk_p(:,jk-2*mk+1)-2*ep^2*xe(:,1).*T.Tk(:,jk-2*mk+1);
            rn0 = find(xe(:,1)>0);
            if p(k)==0
                T.f_prime_over_r(:,k) = (-1)^(jk/2)*(-2*ep^2-jk^2)*ones(n,1);
                T.f_prime_over_r(rn0,k) = -2*ep^2*T.Tk(rn0,jk-2*mk+1)+...
                    T.Tk_p(rn0,jk-2*mk+1)./xe(rn0,1);
                T.f_bis(:,k) = (4*ep^4*T.Pk(:,3) - 2*ep^2).*T.Tk(:,jk-2*mk+1) - ...
                    4*ep^2*T.Pk(:,2).*T.Tk_p(:,jk-2*mk+1) + T.Tk_b(:,jk-2*mk+1);
            else
                T.f_over_r(:,k) = (-1)^((jk-1)/2)*jk*ones(n,1);
                T.f_over_r(rn0,k) = T.Tk(rn0,jk-2*mk+1)./xe(rn0,1);
                T.f_prime_over_r_minus_f_over_r2(rn0,k) = T.Tk_p(rn0,jk-2*mk+1)./xe(rn0,1)-...
                    (2*ep^2+1./T.Pk(rn0,3)).*T.Tk(rn0,jk-2*mk+1);
                T.f_bis(:,k) = (4*ep^4*T.Pk(:,3) - 2*ep^2).*T.Tk(:,jk-2*mk+1) - ...
                    4*ep^2*xe(:,1).*T.Tk_p(:,jk-2*mk+1) + T.Tk_b(:,jk-2*mk+1);
            end
        end
    end
end