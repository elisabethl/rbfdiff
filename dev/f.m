function y = f(op,xk)
x=xk(:,1);
y=xk(:,2);

if (length(op)==1)
  if (op=='0')
    y = sin(x-2*y);
  elseif (op=='x')
    y = cos(x-2*y);
  elseif (op=='y')
    y = -2*cos(x-2*y);
  elseif (op=='L')
    y = -5*sin(x-2*y);
  end
elseif(length(op)==2)
  if (op=='xx')
    y = -sin(x-2*y);
  elseif (op=='xy')
    y =  2*sin(x-2*y);
  elseif (op=='yy')
    y = -4*sin(x-2*y);
  end  
end
