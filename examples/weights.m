function [f s] = weights(psi,ndiff,ptch)

dim = size(ptch.C,2);
%
% Make sure we compute only required derivatives (all of lower order due to quotient rule)
%
ndiffTot = 0;
op = [];
opDim = [];
ndiffVec = [];
while ndiffTot<ndiff
    [opTemp,opDimTemp,~]=getOpDirect(ndiffTot,dim);
    op = [op,opTemp];
    opDim = [opDim, opDimTemp];
    ndiffVec = [ndiffVec, ndiffTot*ones(1,length(opTemp))];
    ndiffTot = ndiffTot + 1;
end
[opTemp,opDimTemp,~]=getOpDirect(ndiff,dim);
op = [op,opTemp];
opDim = [opDim, opDimTemp];
ndiffVec = [ndiffVec, ndiffTot*ones(1,length(opTemp))];

s = zeros(length(unique(vertcat(ptch.xe(:).globalId))),1);
for p = 1:length(ptch.R)
    xLoc = ptch.xe(p).nodes;
    xLoc = (xLoc - ptch.C(p,:))./ptch.R(p); 

    r = xcdist(zeros(1,dim),xLoc,1);
    fTemp = repmat({zeros(size(r))},1,length(op));
    for k=1:length(op)
        fTemp{k} = ((1/ptch.R(p)).^ceil(ndiffVec(k))).*RBFmat(psi,1,r,op{k},opDim{k});
    end
    s(ptch.xe(p).globalId,1) = s(ptch.xe(p).globalId,1) + fTemp{1}'; % Sum of generating function on all patches
    %
    % Organize the output
    %
    if (ndiff==0)
        f{p}.f = fTemp{1};
    elseif (ndiff==1)
        f{p}.f = fTemp{1};
        f{p}.grad = {fTemp{2:end}};
    elseif (ndiff==1.5)
        f{p}.f = fTemp{1};
        f{p}.grad = {fTemp{2:end-1}};
        f{p}.L = fTemp{end};
    elseif (ndiff==2)
        f{p}.f = fTemp{1};
        f{p}.grad = {fTemp{2:dim+1}};
        for k = dim+2:length(op)
            if length(opDim{k}) == 1 
                f{p}.hess{opDim{k},opDim{k}} = fTemp{k};
            else
                f{p}.hess{opDim{k}(1),opDim{k}(2)} = fTemp{k};
                f{p}.hess{opDim{k}(2),opDim{k}(1)} = fTemp{k};
            end
        end
    end
end
%
% Now compute weights using Shepard's method
%
for p = 1:length(ptch.R)
    sLoc = s(ptch.xe(p).globalId,1);
    if (ndiff==0)
        w{p}.f = f{p}.f./sLoc;
    elseif (ndiff==1)
        w{p}.f = fTemp{1}./sLoc;
        w{p}.grad = {fTemp{2:end}};
    elseif (ndiff==1.5)
        w{p}.f = fTemp{1}./sLoc;
        w{p}.grad = {fTemp{2:end-1}};
        w{p}.L = fTemp{end};
    elseif (ndiff==2)
        w{p}.f = fTemp{1}./sLoc;
        w{p}.grad = {fTemp{2:dim+1}};
        for k = dim+2:length(op)
            if length(opDim{k}) == 1 
                w{p}.hess{opDim{k},opDim{k}} = fTemp{k};
            else
                w{p}.hess{opDim{k}(1),opDim{k}(2)} = fTemp{k};
                w{p}.hess{opDim{k}(2),opDim{k}(1)} = fTemp{k};
            end
        end
    end
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

