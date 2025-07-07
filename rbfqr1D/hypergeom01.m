function y=hypergeom01(b,x)
%
% In our case, I will only implement the 0F1 case, which we have.
%
  if (length(b)~=1)
    error('Wrong number of arguments in hypergeom for 0F1')
  end
  %
  % The first coefficient is 1 always
  % 
  alpha=1;
  [ma,pos]=max(abs(x)); % Could possibly be complex
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
    alpha = alpha/(b(1)+n);
    n = n+1;
    %
    % Power and factorial can also be done recursively.
    %
    v = (1/n)*v.*x;
    y = y + alpha*v;
    test = alpha*v(pos);
  end
  
    