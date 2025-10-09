function D=EvalD_3D(ep,p1,p2,j,m,s)
%--- Compute the scaling effect of D_1^{-1} and D_2
%
% Assuming that we get p1 and p2 as vectors.
%
D = zeros(length(p1),length(p2));
%
% Precompute all positive powers of ep that are present.
%
pmax = max(j(p2));
epp = zeros(pmax+1,1);
ep2 = ep*ep;
epp(1)=1;
for p=1:pmax
    epp(p+1) = ep2*epp(p);
end
%
% Precompute powers of 2^4, both negative and postive
%
mmax = max(m(p2));
% mmax = max(j(p2)); % Dirty fix, not sure what to do. Problem relsc!
mmin = max(m(p1)); % Largest negative
twop(0+mmin+1)=1;  % 2^0
for p=1:mmax
    twop(p+mmin+1)=16*twop(p+mmin);
end
for p=1:mmin
    twop(-p+mmin+1)=0.0625*twop(-p+mmin+2);
end
twop2 = [0.5 1 2];
%
% The column values stay the same for each row
%
powj=j(p2);
powm=m(p2);
pows=s(p2);
for k=1:length(p1)
    %
    % Construct powers for one row at a time. Must take away the
    % negative powers
    %
    pow = powj-j(p1(k));
    pos = find(pow>=0);
    D(k,pos) = epp(pow(pos)+1);
    
    pow = powm-m(p1(k));
    D(k,:) = D(k,:).*twop(pow+mmin+1);
    %
    % In the 3D case, we also have a factor 2^(p2-p1), which is either
    % 2,1 or 0.5
    pow = pows-s(p1(k));
    D(k,:) = D(k,:).*twop2(pow+2);
end

f1 = (j(p2)+2*m(p2)+s(p2))/2;
f1t = j(p2)+2*m(p2)+s(p2)+1;
f2 = (j(p2)-2*m(p2)-s(p2))/2;
f3 = (j(p1)+2*m(p1)+s(p1))/2;
f3t = j(p1)+2*m(p1)+s(p1)+1;
f4 = (j(p1)-2*m(p1)-s(p1))/2;
fmax = max(max(f1t),max(f3t));
fp = cumprod([1 1:fmax]); % Because 0!=1
for k=1:length(p1)
    numer = ones(1,size(D,2));
    denom = numer;
    pos = find(f3(k) < f1);
    denom(pos) = (1/fp(f3t(k)+1))*fp(f1t(pos)+1);
    numer(pos) = fp(f1(pos)+1)./fp(f3(k)+1);
    pos = find(f4(k) < f2);
    denom(pos) = (1/fp(f4(k)+1))*denom(pos).*fp(f2(pos)+1);
    pos = find(f3(k) > f1);
    denom(pos) = fp(f3(k)+1)./fp(f1(pos)+1).*denom(pos);
    numer(pos) = fp(f3t(k)+1)./fp(f1t(pos)+1);
    pos = find(f4(k) > f2);
    numer(pos) = numer(pos).*(fp(f4(k)+1)./fp(f2(pos)+1));
    D(k,:) = D(k,:).*numer./denom;
end
