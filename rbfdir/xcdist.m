function r=xcdist(x,c,all)

% if (nargin==1)
%   all=0;
%   c=x;
% elseif (nargin==2)
%   all=0;
% end
% 
% dim = size(x,2);
% switch dim
%     case 3
%         dx = repmat(x(:,1), 1, size(c,1)) - repmat(c(:,1)', size(x,1), 1);
%         dy = repmat(x(:,2), 1, size(c,1)) - repmat(c(:,2)', size(x,1), 1);
%         dz = repmat(x(:,3), 1, size(c,1)) - repmat(c(:,3)', size(x,1), 1);
%         r = sqrt(dx.^2 + dy.^2 + dz.^2);
%     case 2
%         dx = repmat(x(:,1), 1, size(c,1)) - repmat(c', size(x,1), 1);
%         dy = repmat(x(:,2), 1, size(c,1)) - repmat(c', size(x,1), 1);
%         r = sqrt(dx.^2 + dy.^2);
% end


% 
% Evaluation/collocation points and center points are not always the
% same. We compute r=||xi-cj|| together with componentwise signed differences.
%
if (nargin==1)
  all=0;
  c=x;
elseif (nargin==2)
  all=0;
end
  
[np,nd] = size(x);
nc = size(c,1);
nr = 1 + all*nd;

% r contains r and each ri, i=1...nd
r = zeros(np,nc,nr);
if (np==0 | nc==0)
  return
end  
for d=1:nd
  [pi,pj] = meshgrid(c(:,d),x(:,d));
  r(:,:,1) = r(:,:,1) + (pi-pj).^2;
  if (all)
    r(:,:,d+1) = pj-pi;
  end
end
r(:,:,1) = sqrt(r(:,:,1));






