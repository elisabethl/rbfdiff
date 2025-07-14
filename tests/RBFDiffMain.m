clear
close all
% The idea is to test both RBFInterpMat and RBFDiffMat
phi = {'rbfqr','gs','mq'};
pdeg = [-1 -1 3];
dim = 3;
Ndim = [8,28,56];
N = Ndim(dim); % For 1D, 8 is a good number, 2D: 28, 3D:56
Ne = N*4;
epvec = logspace(-2,0,100);
fnum = 1;
if (fnum==1)
    fun = @(x) ones(size(x,1),1);
    funx = @(x) zeros(size(x,1),1);
    funxx = @(x) zeros(size(x,1),1);
end    
%
% We let the domain be a disc or a sphere centered in the origin + 1*n_x
% The node points are Halton nodes without clustering
%
shift = zeros(1,dim); shift(1) = 1;
R = 1;

dimRat = [1 1.2*4/pi 1.2*8/(4*pi/3)];
Nin = ceil(dimRat(dim)*N);
xc = 2*(halton(Nin,dim)-0.5);
r2 = sum(xc.^2,2);
pos = find(r2<=R);
pos = pos(1:N);
xc = xc(pos,:) + shift; 
%
% We do the same for the evaluation points
%
Nin = ceil(dimRat(dim)*Ne);
xe = 2*(halton(Nin,dim)-0.5);
r2 = sum(xe.^2,2);
pos = find(r2<=R);
pos = pos(1:Ne);
xe = xe(pos,:) + shift;

uc = fun(xc);
ue = fun(xe);
uxe = funx(xe);
uxxe = funxx(xe);

for j=1:length(epvec)
    ep = epvec(j);
    for k=1:length(phi)
        %
        % Initialize the approximation
        %
        Psi = RBFInterpMat(phi{k},pdeg(k),ep,xc,shift,R);
        %
        % Compute the differentiation matrices
        %
        [B,Te] = RBFDiffMat(0,Psi,xe);
        [Bx,Te] = RBFDiffMat(1,Psi,Te);
        [BL,Te] = RBFDiffMat(1.5,Psi,Te);
        %
        % Apply them to the right hand side
        %
        u = B*uc;
        Lu = BL*uc;
        %
        % Compute errors
        %
        erru(j,k) = max(abs(u-ue));
        errL(j,k) = max(abs(Lu-uxxe));
    end
end

figure
loglog(epvec,erru)
legend('RBF-QR','GS','MQ+p3')
