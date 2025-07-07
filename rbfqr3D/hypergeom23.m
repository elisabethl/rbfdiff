function y=hypergeom23(a,b,x)
%
% Here, I will only implement the 2F3 case.
%
if (length(a)~=2 || length(b)~=3)
    error('Wrong number of arguments in hypergeom for 2F3')
end
%
% The first coefficient is 1 always
%
alpha=1;
[~,pos]=max(abs(x)); % Could possibly be complex
v = ones(size(x));
y = v;
n=0;
test=1;
mp = eps; % machine precision
while (abs(test) > mp)
    %
    % The hypergeometric coefficient is computed from the arguments.
    % http://en.wikipedia.org/wiki/Hypergeometric_function#The_series_pFq
    %
    alpha = alpha*(a(1)+n)/(b(1)+n)*(a(2)+n)/(b(2)+n)/(b(3)+n);
    n = n+1;
    %
    % Power and factorial can also be done recursively.
    %
    v = (1/n)*v.*x;
    y = y + alpha*v;
    test = alpha*v(pos);
end

