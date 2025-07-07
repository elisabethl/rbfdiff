clear x T R
syms x
jmax = 31; % The highest order we go to
%
% First compute all Chebyshev polynomials up to jmax as well as the difference
%
T{1}=1;
T{2}=x;
for j=2:jmax
    T{j+1} = simplify(2*x*T{j}-T{j-1});
end

for j=1:2:jmax
    R{(j+1)/2} = simplify(T{j+1}-x*diff(T{j+1},x));
end    
%
% Note that high order terms have large coeffs. Perhaps the linear does not hold for long. 
%
