clear
close all
setPaths
%
% Testing the RBFDiffmat implementation for different types of vasis functions
%
dim = 2; % For 1D: 8, 2D: 28, 3D: 56
R = 1;
prob = 4;
switch prob
  case 1
    phi = {'mq','gs','iq'};
    pdeg = [-1 -1 -1]; % One value per phi
    epmin = 2e-1;
    epmax = 5;
    epvec = logspace(log10(epmin),log10(epmax),25);
  case 2
    phi = {'rbfqr'};
    pdeg = [-1]; % One value per phi
    epmin = 1e-2;
    epmax = 2;
    epvec = logspace(log10(epmin),log10(epmax),25);
  case 3
    phi = {'w2'};
    pdeg = [-1]; % One value per phi
    epvec = linspace(0.5,2,10); % The radius of the support
  case 4
    phi = {'phs'};
    pdeg = [3]; % One value per phi
    epvec = [3 5 7 9]; % The order of the PHS
end
lt = {'-','--'};
col = colormap(lines);               
col = col([1 2 4 5 6 7 3],:);

if dim == 1
    Nvec = [8 16];
elseif dim == 2
    Nvec = [28 55];
else
    Nvec = [56 84]; 
end
oversamp = 4; % Oversampling
fnum = 2; % 1 = constant, 2 = sin(pi/2*x*y*z)

syms x y z
X = [x;y;z];

if (fnum==1)
    fun = @(x) ones(size(x,1),1);
    lapFun = @(x) zeros(size(x,1),1);
    for i = 1:dim
        gradFun{i} = @(x) zeros(size(x,1),1);
        for j = 1:dim
            hessFun{i,j} = @(x) zeros(size(x,1),1);
        end
    end
elseif (fnum == 2)
    fun = sin(prod(X(1:dim))*pi/2);
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

nfig = 2 + dim + dim^2;
for k=1:nfig
    fig(k) = figure;
end

for p=1:length(Nvec)
    N = Nvec(p);
    Ne = N*oversamp;
%
% We place halton points on a unit sphere in dim dimensions
%
dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
Nin = ceil(dimRat(dim)*N);
NinE = ceil(dimRat(dim)*Ne);

xc = 2*(halton(Nin,dim)-0.5);
r2 = sqrt(sum(xc.^2,2));
pos = find(r2<=1);
pos = pos(1:N);
xc = xc(pos,:); 
%
% We do the same for the evaluation points
%
xe = 2*(halton(NinE,dim)-0.5);
r2 = sqrt(sum(xe.^2,2));
pos = find(r2<=1);
pos = pos(1:Ne);
xe = xe(pos,:);

if (fnum==1)
    uc = fun(xc);
    ue = fun(xe);
    uL = lapFun(xe);
    for i = 1:dim
        gradU{i} = gradFun{i}(xe);
        for j = 1:dim
            hessU{i,j} = hessFun{i,j}(xe);
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

C = zeros(1,dim); R = 1; % For the unit sphere
for j=1:length(epvec)
    ep = epvec(j);
    for k=1:length(phi)
      phi{k}
        %
        % Initialize the approximation
        %
      Psi = RBFInterpMat(phi{k},pdeg(k),ep,xc,C,R);
        %
        % Compute all differentiation matrices
        %
        [B,~] = RBFDiffMat(0,Psi,xe);
        [BL,~] = RBFDiffMat(1.5,Psi,xe);
        [gradB,~] = RBFDiffMat(1,Psi,xe);
        [hessB,~] = RBFDiffMat(2,Psi,xe);
        erru(j,k) = max(abs(B*uc-ue))/(max(abs(ue))+(max(abs(ue))==0));
        errL(j,k) = max(abs(BL*uc-uL))/(max(abs(uL))+(max(abs(uL))==0));
        for i = 1:dim
            % Compute error in first derivatives
            errGradU{i}(j,k) = max(abs(gradB{i}*uc - gradU{i}))/(max(abs(gradU{i}))+(max(abs(gradU{i}))==0)); 
            for l = 1:dim
                % Compute error in second derivatives
                errHessU{i,l}(j,k) = max(abs(hessB{i,l}*uc - hessU{i,l}))/(max(abs(hessU{i,l}))+(max(abs(hessU{i,l}))==0)); 
            end
        end
    end
end
%
% Plot errors in all derivatives
%
figure(fig(1))
H = loglog(epvec,erru); hold on
fixFig(H,lt(p),col,0)
legend(phi)


figure(fig(2))
H = loglog(epvec,errL); hold on
fixFig(H,lt(p),col,[2 2])
legend(phi)


for i = 1:dim
    figure(fig(2+i))
    H = loglog(epvec,errGradU{i}); hold on
    fixFig(H,lt(p),col,i)
    legend(phi)
end

for i = 1:dim
    for j = 1:dim
        figure(fig(2+dim+((i-1)*dim+j)))
        H = loglog(epvec,errHessU{i,j}); hold on
        fixFig(H,lt(p),col,[i j])
        legend(phi)
    end
end
end
for k=1:length(fig)
    figure(fig(k))
    ax = axis;
    ax(4) = 50; % Cap the y-axis
    axis(ax)
    for j=1:length(Nvec)
        for i=1:length(phi)
            lstr{(j-1)*length(phi)+i} = strcat(phi{i},', N=', num2str(Nvec(j)));
        end    
        H = legend(lstr);
        set(H,'NumColumns',2,'Location','NorthWest')
    end
end

function fixFig(H,lt,col,deriv)
    for k=1:length(H)
        set(H(k),'LineStyle',lt,'Color',col(k,:),'LineWidth',2)  
    end
    set(gca,'FontSize',18)
    set(gcf,'Color','w')
    xlabel('$\varepsilon$','interpreter','latex')
    if deriv==0
        str = 'Evaluation error';
    elseif all(deriv==2)
        str = 'Laplacian error';    
    elseif length(deriv)==1
        str = strcat('Error in $\frac{\partial\phi}{\partial x_',num2str(deriv),'}$');
    elseif deriv(1)==deriv(2)
        str = strcat('Error in $\frac{\partial^2\phi}{\partial x_',num2str(deriv(1)),'^2}$');  
    else    
        str = strcat(' Error in $\frac{\partial^2\phi}{\partial x_',num2str(deriv(1)),'\partial x_',num2str(deriv(2)),'}$');
    end    
    title(str,'interpreter','latex')
end

