clear
close all
setPaths
%
% Testing the RBFmat implementation
%
phi = {'r3','rbfqr','mq','gs'};
pdeg = [3 -1 3 2];
dim = 3; % For 1D: 8, 2D: 28, 3D: 56
if dim == 1
    N = 8;
elseif dim == 2
    N = 28;
else
    N = 56; 
end
Ne = N*4; % Oversampling
epvec = logspace(-2,0,100);
fnum = 1;

syms x y z
X = [x;y;z];

if (fnum==1)
    fun = prod(X(1:dim)).^3;
    lapFun = 0;
    for i = 1:dim
        gradFun{i} = matlabFunction(diff(fun,X(i)));
        lapFun = lapFun + diff(fun,X(i),X(i));
        for j = 1:dim
            hessFun{i,j} = matlabFunction(diff(fun,X(i),X(j)));
        end
    end
    fun = matlabFunction(fun);
    lapFun = matlabFunction(lapFun);
elseif (fnum == 2)
    fun = sin(prod(X(1:dim))*pi);
    lapFun = 0;
    for i = 1:dim
        gradFun{i} = matlabFunction(diff(fun,X(i)));
        lapFun = lapFun + diff(fun,X(i),X(i));
        for j = 1:dim
            hessFun{i,j} = matlabFunction(diff(fun,X(i),X(j)));
        end
    end
    fun = matlabFunction(fun);
    lapFun = matlabFunction(lapFun);
end    
%
% We place halton points in a circle centred at "shift" and stretched with R
%
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
Nin = ceil(dimRat(dim)*N);
NinE = ceil(dimRat(dim)*Ne);
shift = [-1,3.2,-4]; % Point centre in 3D case
shift = shift(1:dim);
R = 2;
xc = 2*R*(halton(Nin,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=R);
pos = pos(1:N);
xc = xc(pos,:) + shift; 
%
% We do the same for the evaluation points
%
xe = 2*R*(halton(NinE,dim)-0.5);
r2 = sqrt(sum(xe.^2,2));
pos = find(r2<=R);
pos = pos(1:Ne);
xe = xe(pos,:) + shift;

if (fnum==1)
    cP = num2cell(xc,1);
    uc = fun(cP{:});
    cP = num2cell(xe,1);
    ue = fun(cP{:});
    uL = lapFun(cP{:});
    for i = 1:dim
        gradU{i} = gradFun{i}(cP{:});
        for j = 1:dim
            hessU{i,j} = hessFun{i,j}(cP{:});
        end
    end
elseif (fnum == 2)
    cP = num2cell(xc,1);
    uc = fun(cP{:});
    cP = num2cell(xe,1);
    ue = fun(cP{:});
    uL = lapFun(cP{:});
    for i = 1:dim
        gradU{i} = gradFun{i}(cP{:});
        for j = 1:dim
            hessU{i,j} = hessFun{i,j}(cP{:});
        end
    end
end  

for j=1:length(epvec)
    ep = epvec(j);
    for k=1:length(phi)
        %
        % Initialize the approximation
        %
        Psi = RBFInterpMat(phi{k},pdeg(k),ep,xc,shift,R);
        %
        % Compute all differentiation matrices
        %
        [B,~] = RBFDiffMat(0,Psi,xe);
        [BL,~] = RBFDiffMat(1.5,Psi,xe);
        [gradB,~] = RBFDiffMat(1,Psi,xe);
        [hessB,~] = RBFDiffMat(2,Psi,xe);
        erru(j,k) = max(abs(B*uc-ue))/max(abs(ue));
        errL(j,k) = max(abs(BL*uc-uL))/max(abs(uL));
        for i = 1:dim
            % Compute error in first derivatives
            errGradU{i}(j,k) = max(abs(gradB{i}*uc - gradU{i}))/max(abs(gradU{i})); 
            for l = 1:dim
                % Compute error in second derivatives
                errHessU{i,l}(j,k) = max(abs(hessB{i,l}*uc - hessU{i,l}))/max(abs(hessU{i,l})); 
            end
        end
    end
end

figure
loglog(epvec,erru)
legend(phi)
title("Evaluation error")

figure
loglog(epvec,errL)
legend(phi)
title("Laplace error")

for i = 1:dim
    figure()
    loglog(epvec,errGradU{i});
    legend(phi)
    title(['dx',num2str(i),' error']) 
end

for i = 1:dim
    for j = 1:dim
        figure()
        loglog(epvec,errHessU{i,j});
        legend(phi)
        title(['dx',num2str(i),num2str(j),' error']) 
    end
end
