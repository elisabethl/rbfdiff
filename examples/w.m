function weights = w(psi,ndiff,ptch,x,ptPtchLst)
dim = size(ptch.C,2);
compTrue = 0;
if (ndiff > 0 & ndiff~=1.5) % Evaluation or Laplacian
    compTrue = 1; % Meaning that we need component distances as well
end    
s = zeros(size(x,1),1);
for p = 1:length(ptch.R)
    idLoc = find(ptPtchLst==p);
    xLoc = x(idLoc,:);
    xLoc = (xLoc - ptch.C(p,:))./ptch.R; 

    r = xcdist(zeros(1,dim),xLoc,compTrue);
    [op,opDim,derVec]=getOpDirect(ndiff,dim);
    s = RBFmat(psi,1,r,'0',1);
    
    for k=1:length(op)
        f{p} = ((1/ptch.R(p)).^ceil(sum(ndiff))).*RBFmat(psi,1,r,op{k},opDim{k});
    end
    sTot(idLoc,1) = sTot(idLoc,1) + s;
end
%
% Now compute weights using Shepard's method
%
for p = 1:length(ptch.R)
    idLoc = find(ptPtchLst==p);
    xLoc = x(idLoc,:);
    xLoc = (xLoc - ptch.C(p,:))./ptch.R; 

    r = xcdist(zeros(1,dim),xLoc,compTrue);
    [op,opDim,derVec]=getOpDirect(ndiff,dim);
    
    for k=1:length(op)
        f = RBFmat(psi,1,r,op{k},opDim{k});
    end
    s(idLoc,1) = s(idLoc,1) + f;
end

function [op,opDim,derVec]=getOpDirect(ndiff,dim)
I = eye(dim);
if (ndiff==0)
    op{1} = '0';
    opDim{1} = 1;
    derVec{1} = zeros(1,dim);
elseif (ndiff==1)
    for k=1:dim
        op{k} = '1';
        opDim{k} = k;
        derVec{k} = I(k,:);
    end
elseif (ndiff==1.5)
    op{1} = 'L';
    opDim{1} = dim;
    derVec{1} = 2*I;
elseif (ndiff==2)
    pos = 1;
    for k=1:dim
        op{pos} = '2';
        opDim{pos} = k;
        derVec{pos} = 2*I(k,:);
        pos = pos + 1;
        for j=k+1:dim
            op{pos} = 'm2';
            opDim{pos} = [j k];
            derVec{pos} = I(k,:) + I(j,:);
            pos = pos + 1;
        end    
    end
elseif (ndiff > 2)
    if mod(ndiff,2)==0 % Positive values for higher order Laplacian
                       % We could also encode them with half numbers
        op{1} = strcat('L',num2str(ndiff/2));
        opDim{k} = 1;
        %
        % Add code here
        %
    end    
end

