% Computing (T_j-r*T_j^\prime)/r^2 for odd j
function R = Rjm0(T,dT,r)
%
% An expansion is used for small x and direct evaluation for large x.
% The breakpoint is determined such that the rounding error  for direct
% evaluation is around 5e-13 (2*eps*jmax./breakpt) for jmax = 51, which is
% in general larger than what makes sense to use.
% For double precision, nmax terms are then enough for the expansion. Should this code be used in higher precision, nmax can be adjusted or set to jend
% (must be an odd number).
%
breakpt = 0.04;
nmax = 19;
r = r(:); % Ensure column vector
r2 = r.^2;
jmax = size(T,2)-1;
jvec = 3:2:jmax; % Only odd indices from 3 on have a non-zero value
jhalf = 0.5*(jvec+1);
jend = jvec(end); % The last odd index
%
% Divide r into small and large values
%
Rl = zeros(0,length(jvec));
Rs = zeros(0,length(jvec));
isSmall = (abs(r)<=breakpt);
psmall = find(isSmall);
plarge = find(~isSmall);
%
% Use direct evaluation for large values
%
if ~isempty(plarge)
    Rl = (T(plarge,jvec+1)-r(plarge).*dT(plarge,jvec+1)).*(1./r2(plarge));
end    
%
% Use a recursive Horner-like scheme for small values
%
x = r(psmall); x2 = r2(psmall);
nterm = min(nmax,jend); % Must be an odd number
%
if ~isempty(psmall)
    Rs = (nterm-1)/factorial(nterm); 
    for k=nterm-2:-2:1
        Rs = -Rs.*(x.^2*(jvec.^2-k^2)) + (k-1)/factorial(k);  
    end
    Rs = Rs.*(1./x*(jvec.*(-1).^jhalf)); % Alternating signs
end
%
% Assign to desired global indices
%
R = zeros(length(r),jmax+1);
R(psmall,jvec+1) = Rs;
R(plarge,jvec+1) = Rl;
