function [phi]=si(epsil,r,nprime,dim)

if nargin<=3
  dim=1;
end

phi=zeros(size(r,1),size(r,2));

tmp = epsil*sq(r(:,:,1));

if nprime(1)=='0'

  phi = sinc(tmp/pi);

elseif nprime(1)=='L' & length(nprime)==1

nd = size(r,3)-1;

mask=sq(r(:,:,1))==0; 

phi = ((nd-3)*(tmp.*cos(tmp)-sin(tmp))-tmp.^2.*sin(tmp))./ ...
      (epsil*sq(r(:,:,1)).^3+mask);
phi=phi.*(mask==0);

phi = phi-nd/3*epsil^2*mask;

elseif nprime(1:2)=='L2'

nd = size(r,3)-1;

mask=sq(r(:,:,1))==0; 
phi = ( (-3*nd^2+24*nd-45)*(tmp.*cos(tmp)-sin(tmp)) + ...
        (-nd^2+10*nd-21)*tmp.^2.*sin(tmp) + ...
      (-2*nd+6)*tmp.^3.*cos(tmp) + tmp.^4.*sin(tmp) )./(epsil*sq(r(:,:,1)).^5+mask);
phi=phi.*(mask==0);
phi = phi + nd/15*(nd+2)*epsil^4*mask;

else
  error('Error in input argument nprime to function si')
end 

function r=sq(r)
 r=squeeze(r);
