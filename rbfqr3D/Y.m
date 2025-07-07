%
% For a description of the quantities that are computed, see the article
% [] and the handout on second order derivatives for RBF-QR in 3D.
% The last two derivatives include combinations t remove singularities.
%
function [YY,Yth,Yfi_over_s,Yfi,Ythth,Ythfi,Yfifi]=Y(mu,th,fi,YY)
%%% [YY,Yth,Yfi_over_s] = Y(mu,th,fi,YY)
%
% Given a max SPH level mu and column vectors with th and fi
% locations on the unit sphere, function Y returns a
% (length(th))x(mu+1)x(mu+1) array with all SPH values up to level
% mu at all the spatial locations.
% In the top (1,:,:)-plane of YY, the Y(mu,nu)-entries are located
% as follows:
%          [(0, 0) (1, 1) (2, 2) (3, 3) ....]
%          [(1,-1) (1, 0) (2, 1) (3, 2) ....]
%          [(2,-2) (2,-1) (2, 0) (3, 1) ....]
%          [(3,-3) (3,-2) (3,-1) (3, 0) ....]
%          [ ....   ....   ....   ....  ....]
%
if nargin>3 && ~isempty(YY)
    YYold = YY;
    mmin = size(YY,2)-1;
    YY = zeros(length(th),mu+1,mu+1);
    YY(:,1:mmin+1,1:mmin+1) = YYold;
else
    YY = zeros(length(th),mu+1,mu+1);
    mmin = 0;
end

if nargout > 1
    Yth = zeros(length(th),mu+1,mu+1);
    Yfi_over_s = zeros(length(th),mu+1,mu+1);
end

if nargout > 3
    Yfi   = zeros(length(th),mu+1,mu+1);
    Ythth = zeros(length(th),mu+1,mu+1);
    Ythfi = zeros(length(th),mu+1,mu+1);
    Yfifi = zeros(length(th),mu+1,mu+1);
end
% Do case (0,0) separately (Matlab's legendre routine a bit strange
% in this special case).
c = 1/sqrt(2*pi);
YY(:,1,1) = c*sqrt(0.5);
z = cos(th);
%
% We need the inverse of sine theta for some second derivatives.
% The limit values for theta = 0, pi are zero.
%
        
[t1,t2]=meshgrid(1:mu,fi); t = t1.* t2;
cosfi = cos(t); sinfi = sin(t);
for m=mmin+1:mu
    %   p = legendre(m,z,'norm')'; % We need the normalization 'norm'.
    [p,pth,p_over_s,pthth,pfifi,pthfi] = legendre_modified(m,z,'norm');
    p = p.'; pth = pth.'; p_over_s = p_over_s.'; pthth = pthth.';
    pfifi = pfifi.'; pthfi=pthfi.';
    YY(:,m+1,m+1) = c*p(:,1);
    for k=1:m
        YY(:,m+1-k,m+1) = c*p(:,k+1).*cosfi(:,k);
        YY(:,m+1,m+1-k) = c*p(:,k+1).*sinfi(:,k);
    end
    if nargout > 1
        Yth(:,m+1,m+1) = -c*pth(:,1);
        for k=1:m
            Yth(:,m+1-k,m+1) = -c*pth(:,k+1).*cosfi(:,k);
            Yth(:,m+1,m+1-k) = -c*pth(:,k+1).*sinfi(:,k);
            Yfi_over_s(:,m+1-k,m+1) = -c*k*p_over_s(:,k+1).*sinfi(:,k);
            Yfi_over_s(:,m+1,m+1-k) = c*k*p_over_s(:,k+1).*cosfi(:,k);
        end
    end
    if (nargout > 3)
        %
        % Note that fifi and thfi are combinations, not pure derivatives.
        % The mixed second derivative is zero for nu=0;
        %
        Ythth(:,m+1,m+1) = c*pthth(:,1);
        Yfifi(:,m+1,m+1) = c*pfifi(:,1);
        for k=1:m
          Yfi(:,m+1-k,m+1) = -c*k*p(:,k+1).*sinfi(:,k);
          Yfi(:,m+1,m+1-k) =  c*k*p(:,k+1).*cosfi(:,k);
          
          Ythth(:,m+1-k,m+1) = c*pthth(:,k+1).*cosfi(:,k);
          Ythth(:,m+1,m+1-k) = c*pthth(:,k+1).*sinfi(:,k);
          
          Yfifi(:,m+1-k,m+1) = c*pfifi(:,k+1).*cosfi(:,k);
          Yfifi(:,m+1,m+1-k) = c*pfifi(:,k+1).*sinfi(:,k); 

          Ythfi(:,m+1-k,m+1) = c*pthfi(:,k+1).*sinfi(:,k);
          Ythfi(:,m+1,m+1-k) = -c*pthfi(:,k+1).*cosfi(:,k);
        end  
    end  
end