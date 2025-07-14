%
% The purpose of this script is to test the RBFDiffMat computes correct
% derivatives for sensible parameter combinations. The tester can choose
% basis function and dimension and will be presented with convergence plots.
%
% TODO: Calculations, and perhaps give x and xe as inparameters. Then we can try special regions. Perhaps do one for 1D and one for 2D etc, but these are in a sense the mains that we have. We want to be able to have '1'
function [epvec,relerr] = testRBFDiffMat(phi,dim,fnum)
close all

x = sym('x');
y = sym('y');
z = sym('z');

f = cos(pi/3*x)*sin(pi*y)*cos(0.5*pi*z); z0 = 0; y0 = 0.5; % To be one
f = matlabFunction(f);
fx = matlabFunction(diff(f,x));
fxx = matlabFunction(diff(fx,x));
fxy = matlabFunction(diff(fx,y));
fxz = matlabFunction(diff(fx,z));
fy = matlabFunction(diff(f,y));
fyy = matlabFunction(diff(fy,y));
fyz = matlabFunction(diff(fy,z));
fz = matlabFunction(diff(f,z));
fzz = matlabFunction(diff(fz,z));
fLap{1} = fxx; % Cannot plus here
fLap{2} = @(x,y,z) fxx(x,y,z) + fyy(x,y,z);
fLap{3} = @(x,y,z) fxx(x,y,z) + fyy(x,y,z) + fzz(x,y,z);

fgrad = {fx,fy,fz};
fHess= {fxx,fxy,fxz;0,fyy,fyz;0,0,fzz};

% All derivatives that are tested in the first round
igrad = [2 5 7];
iHess = [3 4 9 ]
flist = {f,fx,fxx,fxy,fy,fyy,fz,fzz,fxz,fyz,fLap2,fLap3};
fnum=[3 6 10];
derVec = [0 0 0; %f
          1 0 0; %fx
          2 0 0; %fxx
          1 1 0; %fxy
          0 1 0; %fy
          0 2 0; %fyy
          0 0 1; %fz
	  0 0 2, %fzz
	  1 0 1; %fxz
          0 1 1]; %fyz
fstr={'f','fx','fxx','fxy','fy','fyy','fz','fzz','fxz','fyz','L'};
if (dim==2)
    fstr={'f','fx','fxx','fxy','fy','fyy','L'};
end        
% With or without a polynomial term
if (strcmp(phi,'rbfqr'))
    % We do not add polynomials in the RBF-QR case, since the basis functions
    % themselves are close to polynomials in the limit
    pdegv = -1; 
else    
    pdegv = [-1 3];
end
    

nodeGen = 'halt';

% To make convergence cuves, we select a few different N and a range of
% shape parameters
if strcmp(phi,'rbfqr')
    epvec = logspace(-2,log10(2),25);
else
    epvec = logspace(log10(0.25),log10(5),25);
end    
if (dim==1)
    nVec =[4 7 10 15];
elseif (dim==2)
    nVec = [28 36 45 55];
elseif (dim==3)
    nVec = [35 56 84 120];
end    

if (dim==1)
  argstr = '(xc(:,1),y0,z0);';
elseif (dim==2)
  argstr = '(xc(:,1),xc(:,2),z0);';
else
  argstr = '(xc(:,1),xc(:,2),xc(:,3));';
end

for k=1:length(nVec)
    nLoc = nVec(k)
    % We only test for Halton nodes in a cube
    xc = halton(nLoc,dim);
    xe = halton(2*nLoc,dim);
    fun = matlabFunction(f);
    fstr = strcat('fc=fun',argstr);
    eval(fstr)

    for j=1:length(epvec)
        ep = epvec(j);
        for ell=1:length(pdegv)
            [Psi] = RBFInterpMat(phi,pdegv(ell),ep,xc,xe);
            [E,Te] = RBFDiffMat(0,Psi,xe);
	    [G,Te] = RBFDiffMat(1,Psi,Te);
	    [H,Te] = RBFDiffMat(2,Psi,Te);
	    [L,Te] = RBFDiffMat(1.5,Psi,Te);
	    ind = 1;
	    fderNum(ind,:) = B*fc;
	    fstr = strcat('fderEx(ind,:) = fun',argstr);
	    eval(fstr);
	    ind = ind+1;
            for m=1:dim
	      fderNum(ind,:) = G{m}*fc;
	      fstr = strcat('fderEx(ind,:) = fgrad{m}',argstr);
	    end  
	    
            for m=1:1+dim+dim*(dim+1)/2
	    
	        dfun = matlabFunction(flist{m});
	        if (dim==1)
	            fderEx = dfun(xe(:,1),y0,z0);                    
	        elseif (dim==2)
	            fderEx = dfun(xe(:,1),xe(:,2),z0);
	        else
	            fderEx = dfun(xe(:,1),xe(:,2),xe(:,3));
	        end
	        relerr(j,k,m,ell) = max(abs(fderNum-fderEx))/max(abs(fderEx));
            end
            %
            % Also add a test for the Laplacian if dim > 1
            %
            if (dim > 1)
                derLap = zeros(1,dim); derLap(1) = -1;
                [B,Psi] = RBFDiffMat(derLap,Psi,xe);
	        fderNum = B*fc;
	        if (dim==2)
                    dfun = matlabFunction(flist{end-1});
                    fderEx = dfun(xe(:,1),xe(:,2),z0);
	        elseif (dim==3)
                    dfun = matlabFunction(flist{end});
	            fderEx = dfun(xe(:,1),xe(:,2),xe(:,3));
	        end
	        relerr(j,k,fnum(dim)+1,ell) = max(abs(fderNum-fderEx))/max(abs(fderEx));
            end
        end
    end
end

for m=1:fnum(dim)+1
    for ell=1:length(pdegv)
        figure, H = loglog(epvec,relerr(:,:,m,ell),'-');
        set(H,'LineWidth',2)
	title(['phi=' phi ', dim=' num2str(dim) ', ' fstr{m} ', pdeg=' num2str(pdegv(ell))])
        legend(num2str(nVec(:)))
    end
end    




