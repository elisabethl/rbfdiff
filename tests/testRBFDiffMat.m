%
% The purpose of this script is to test that RBFDiffMat computes correct
% derivatives for sensible parameter combinations. The tester can choose
% basis function and dimension and will be presented with convergence plots.
%
% TODO: Calculations, and perhaps give x and xe as inparameters. Then we can try special regions. Perhaps do one for 1D and one for 2D etc, but these are in a sense the mains that we have. We want to be able to have '1'
function [epvec,relerr,relgrad,relhess,rellap] = testRBFDiffMat(phi,dim,fnum)
close all

x = sym('x');
y = sym('y');
z = sym('z');

switch fnum
  case 1
    f = cos(pi/3*x)*sin(pi*y)*cos(0.5*pi*z); z0 = 0; y0 = 0.5; % To be one
  case 2  
    f = (x+y+z)^0; y0=1; z0=1;
  case 3
    f = 1/(16+x^2 + y^2 + z^2); y0 = 0; z0 = 0; 
  otherwise
    disp('The requested function is not implemented, pls adjust fnum')
    return
end
f = matlabFunction(f,'Vars',{'x','y','z'});
fx = matlabFunction(diff(f,x),'Vars',{'x','y','z'});
fxx = matlabFunction(diff(fx,x),'Vars',{'x','y','z'});
fxy = matlabFunction(diff(fx,y),'Vars',{'x','y','z'});
fxz = matlabFunction(diff(fx,z),'Vars',{'x','y','z'});
fy = matlabFunction(diff(f,y),'Vars',{'x','y','z'});
fyy = matlabFunction(diff(fy,y),'Vars',{'x','y','z'});
fyz = matlabFunction(diff(fy,z),'Vars',{'x','y','z'});
fz = matlabFunction(diff(f,z),'Vars',{'x','y','z'});
fzz = matlabFunction(diff(fz,z),'Vars',{'x','y','z'});
fLap{1} = fxx; % Cannot plus here
fLap{2} = @(x,y,z) fxx(x,y,z) + fyy(x,y,z);
fLap{3} = @(x,y,z) fxx(x,y,z) + fyy(x,y,z) + fzz(x,y,z);

fgrad = {fx,fy,fz};
fHess= {fxx,fxy,fxz;
        0,fyy,fyz;
        0,0,fzz};

fnme={'$f$','$f_x$','$f_y$','$f_z$'}
hstr = {'$f_{xx}$'};
if (dim==2)
  hstr =   {'$f_{xx}$','$f_{xy}$','$f_{yy}$'};
else
  hstr =   {'$f_{xx}$','$f_{xy}$','$f_{xz}$','$f_{yy}$','$f_{yz}$','$f_{zz}$'};
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

% To make convergence curves, we select a few different N and a range of
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
  arge = '(xe(:,1),y0,z0);';
elseif (dim==2)
  argstr = '(xc(:,1),xc(:,2),z0);';
  arge = '(xe(:,1),xe(:,2),z0);';
else
  argstr = '(xc(:,1),xc(:,2),xc(:,3));';
  arge = '(xe(:,1),xe(:,2),xe(:,3));';
end

for k=1:length(nVec)
    nLoc = nVec(k)
    % We only test for Halton nodes in a cube
    xc = halton(nLoc,dim);
    xe = halton(2*nLoc,dim);
    % Evaluate the function at the node points
    fstr = strcat('fc=f',argstr);
    eval(fstr)
    fc = fc.*ones(size(xc(:,1))); % Such that the constant becomes a vector

    for j=1:length(epvec)
        ep = epvec(j);
        for ell=1:length(pdegv)
            [Psi] = RBFInterpMat(phi,pdegv(ell),ep,xc,xe);
            [E,Te] = RBFDiffMat(0,Psi,xe);
	    [G,Te] = RBFDiffMat(1,Psi,Te);
	    [H,Te] = RBFDiffMat(2,Psi,Te);
	    [L,Te] = RBFDiffMat(1.5,Psi,Te);
	    %
	    % Evaluate the interpolant and corresponding error
	    %
	    fnum = E*fc;
	    fstr = strcat('fe=f',arge);
	    eval(fstr)
	    relerr(j,k,ell) = max(abs(fnum-fe))/max(abs(fe));
	    %
	    % Compute numerical gradients and their errors
	    % 
            for d=1:dim
	      gnum = G{d}*fc;
	      fstr = strcat('gev=fgrad{',num2str(d),'}',arge);
	      eval(fstr)
	      relgrad(j,k,ell,d) = max(abs(gnum-gev))/max(abs(gev)+all(gev==0));
	    end
            %
	    % Compute Hessian and errors
	    %
	    for d=1:dim
	      for m=d:dim
		hnum = H{d,m}*fc;
		fstr = strcat('he=fHess{d,m}',arge);
		eval(fstr)
		relhess(j,k,ell,d,m) = max(abs(hnum-he))/max(abs(he)+all(he==0));
	      end
	    end
            %
	    % Compute Laplacian and error
	    %
            lnum = L*fc;
	    fstr = strcat('lev=fLap{',num2str(dim),'}',arge);
	    eval(fstr)
	    rellap(j,k,ell) = max(abs(lnum-lev))/max(abs(lev)+all(lev==0)); 
	end
    end
end

for ell=1:length(pdegv)
  figure, Q = loglog(epvec,relerr(:,:,ell),'-');
  funfix(Q,phi,dim,fnme{1},pdegv(ell),nVec)

  for d=1:dim
    figure, Q = loglog(epvec,relgrad(:,:,ell,d),'-');
    funfix(Q,phi,dim,fnme{1+d},pdegv(ell),nVec)
  end  

  p = 0;
  for d=1:dim
    for m=d:dim
      p=p+1;
      figure, G = loglog(epvec,relhess(:,:,ell,d,m),'-');
      funfix(G,phi,dim,hstr{p},pdegv(ell),nVec)
    end
  end
end

function funfix(H,phi,dim,fstr,pdeg,nVec)
  set(H,'LineWidth',2)
  set(gca,'FontSize',18,'LineWidth',1.5)
  title([fstr ', $\varphi=$ ' phi ', $d=' num2str(dim) '$,  $p=' num2str(pdeg) '$'],'interpreter','latex')
  lstrs = [repmat('$n = ',length(nVec),1) num2str(nVec(:))];
  lstrs(:,end+1) = '$';
  legend(lstrs,'interpreter','latex')
  set(gcf,'Color','w')
  xlabel('$\varepsilon$','interpreter','latex')
  ylabel('Relative error','interpreter','latex')
