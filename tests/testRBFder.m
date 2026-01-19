%
% For the RBF phi with shape parameter ep, the RBF and all its available
% derivatives are plotted in dim dimensions.
% The results are compared with finite difference approximations from the
% previous derivative level e.g. d^4 f/dx_1^4 = d/dx_1(d^3 f/dx_1^3).
% For 3D, the relative error between the 2 is instead provided as a full
% visual comparison is difficult.
function phival=testRBFder(phi,ep,dim)
close all  
%
% We let the domain be [-1 1]^dim with the RBF centered at the origin.
% For now, only dimensions up to 3 are supported.
%
xc = zeros(1,dim);
%
% We choose different numbers of points per dimension depending on dim.
%
Ndim = [201, 101, 41];
x = linspace(-1,1,Ndim(dim));
h = x(2)-x(1); % Step length in finite difference approximation
if dim==1
    xe = x(:);
elseif dim==2
    [xx,yy]=meshgrid(x);
    xe = [xx(:) yy(:)];
elseif dim==3
    [xx,yy,zz]=meshgrid(x);
    xe = [xx(:) yy(:) zz(:)];
else
    error('Only dimensions 1--3 are supported')
end
%
% Compute all possible derivatives, with some limitations for phs and w2
%
dermax = 4;
if strcmp(phi,'w2')
    dermax = 2;
elseif strcmp(phi,'phs')
    dermax = ep;
end
%
% Set up the derivatives and the corresponding arguments to phi.m
%
nprime{1} = '0'; darg{1} = 1; dervec(1,:) = zeros(1,dim); str{1} = '$\phi(r)$'; 
for k=1:dermax
    [np,da,dv,st] = getPars(k,dim);
    nprime = [nprime; np];
    darg = [darg; da];
    dervec = [dervec;dv];
    str = [str;st];
end
%
% Compute the distance matrix \|xe-xc\|
% 
all = 1;
r = xcdist(xe,xc,all); % Compute also component distances
%
% Go through all cases and plot the analytical and FD values
%
for k=1:length(nprime)
    figure
    %
    % Compute and plot the values for each derivative
    %
    phival(:,k) = feval(phi,ep,r,nprime{k},darg{k});
    phiplot = phival(:,k);
    if (dim==1)
        plot(xe,phiplot)
    elseif dim==2
        phiplot=reshape(phiplot,size(xx));
        H = mesh(xx,yy,phiplot);
        set(H,'FaceColor','interp')
    else
        phiplot=reshape(phiplot,size(xx));
        H=slice(xx,yy,zz,phiplot,0,0,0);
        set(H,'FaceColor','interp','EdgeColor','none')
        axis equal
    end
    hold on
    %
    % Now let's make a finite difference approximation if that is possible
    %
    meshdims = [1 0 0;
                2 1 0;
                2 1 3]; % Where to diff given a certain diemension and coord
    [maxd,loc] = max(dervec(k,:));
    if (maxd>0)
      targ = dervec(k,:); targ(loc) = maxd-1;  
      [~,row]=ismember(targ,dervec,'rows');
      fdapprox = phival(:,row);
      if (dim>1)
          sz = Ndim(dim)*ones(1,dim);
          fdapprox = reshape(fdapprox,sz);
      end
      
      fdapprox = 1/h*diff(fdapprox,1,meshdims(dim,loc));
      %
      % We want to compare that with the analytical value at the same grid
      % We compute also that using a difference, so not maximal accuracy
      %
      phicmp = 0.5*diff(phiplot,1,meshdims(dim,loc));
      
      if (dim==1)
          plot(xe(1:end-1)+h/2,fdapprox)
      elseif dim==2
          if (loc==1)
              H = mesh(xx(:,1:end-1)+h/2,yy(:,1:end-1),fdapprox);
          else
              H = mesh(xx(1:end-1,:),yy(1:end-1,:)+h/2,fdapprox);
          end
          set(H,'FaceColor','none','EdgeColor',[0.7 0.7 0.7])
      else
          mid = ceil(Ndim(dim)/2);
          fig = figure; % Just to have somewhere to dump the contours first
          if (loc==1)
              phicmp = phicmp + phiplot(:,1:end-1,:);
              [C1,H] = contour(xx(:,1:end-1,mid)+h/2,yy(:,1:end-1,mid), ...
                               fdapprox(:,:,mid));
              [C2,H] = contour(squeeze(zz(mid,1:end-1,:)), ...
                               squeeze(xx(mid,1:end-1,:))+h/2, ...
                               squeeze(fdapprox(mid,:,:)));
              % This one is a little bit off in terms of where the value is
              [C3,H] = contour(squeeze(zz(:,mid-1,:)), ...
                               squeeze(yy(:,mid-1,:)), ...
                               squeeze(fdapprox(:,mid-1,:)));

          elseif (loc==2)
              phicmp = phicmp + phiplot(1:end-1,:,:);
              [C1,H] = contour(xx(1:end-1,:,mid),yy(1:end-1,:,mid)+h/2, ...
                               fdapprox(:,:,mid));
              [C2,H] = contour(squeeze(zz(mid-1,:,:)), ...
                               squeeze(xx(mid-1,:,:)), ...
                               squeeze(fdapprox(mid,:,:)));
              % This one is a little bit off in terms of where the value is
              [C3,H] = contour(squeeze(zz(1:end-1,mid,:)), ...
                               squeeze(yy(1:end-1,mid,:))+h/2, ...
                               squeeze(fdapprox(:,mid,:)));
          else
              phicmp = phicmp + phiplot(:,:,1:end-1);
              [C1,H] = contour(xx(:,:,mid-1),yy(:,:,mid-1), ...
                               fdapprox(:,:,mid-1));
              [C2,H] = contour(squeeze(zz(mid,:,1:end-1))+h/2, ...
                               squeeze(xx(mid,:,1:end-1)), ...
                               squeeze(fdapprox(mid,:,:)));
              % This one is a little bit off in terms of where the value is
              [C3,H] = contour(squeeze(zz(:,mid,1:end-1))+h/2, ...
                               squeeze(yy(:,mid,1:end-1)), ...
                               squeeze(fdapprox(:,mid,:)));
          end
          close(fig) % Temporary figure
                     % Now plot the contours in the 3D-figure on each slice.
          plotContour(C1,[1 2 3])
          plotContour(C2,[3 1 2])
          plotContour(C3,[3 2 1])
          relerr(k,1) = max(max(max(abs(fdapprox-phicmp))))./ ...
                        max(max(max(abs(phicmp))));
      end
    end
    set(gca,'FontSize',18)
    title(str{k},'interpreter','latex')
    xlabel('$x_1$','interpreter','latex')
    if dim>=2
        ylabel('$x_2$','interpreter','latex')
    end
    if k>1 & dim==1
        legend('Exact','FD')
    end    
    if k>1 & dim==2 % Not for the first function value
        legend('Exact','FD','Location','SouthOutside','Orientation','Horizontal')
    end
    if (dim==3)
        zlabel('$x_3$','interpreter','latex')
        colorbar
    end    
end
if (dim==3)
    disp('In the figures, the colored planes show the analytical results')
    disp('and the white lines are finite difference contours.')
    disp(' ')
    disp('The relative maximum error between the analytical and fd results' )
    disp(num2str([(2:length(relerr))', round(relerr(2:end),2,'significant')]))
end    
% Note that L and L2 have not been tested here

function [nprime,darg,derMat,str] = getPars(der,dim)
%
% First collect all combinations for the given number of dimensions
%
if dim==1
    derMat(1,:) = der;
elseif dim==2
    derMat = [(der:-1:0)' (0:der)'];    
else
    i=0;
    for j=der:-1:0
        for k=(der-j):-1:0
            i = i+1;
            derMat(i,:) = [j k der-j-k];
        end
    end    
end
for k=1:size(derMat,1)
    curr = derMat(k,:);
    if max(curr)==der
        nprime{k,1} = num2str(der);
        darg{k,1} = find(curr==der); % Should be only one position
        pstr = '';
        if der>1
            pstr = strcat('^',num2str(der));
        end    
        str{k,1} = strcat('$\frac{\partial', pstr, ...
                          '\phi}{\partial x_',num2str(darg{k,1}),pstr,'}$');
    else
        nprime{k,1} = strcat('m',num2str(der));
        darg{k,1}=[];
        str{k,1} = strcat('$\frac{\partial^',num2str(der),'\phi}{');
        for i=1:dim
            darg{k,1}=[darg{k,1} i*ones(1,curr(i))];
            if curr(i)>0
                if curr(i)==1
                    pstr = '';
                else
                    pstr = strcat('^',num2str(curr(i)));
                end    
                str{k,1} = strcat(str{k,1},'\partial x_',num2str(i),pstr); ...
            end
        end
        str{k,1} = strcat(str{k,1},'}$');
    end    
end

% Let dims be the two that come in as X and Y, and the third, which will be zero
function plotContour(C,dims)
C(end+1,:) = 0; % We add a row of zeros
pos = 1;
while pos <= size(C,2)
    % Read the number of points
    Ni=C(2,pos);
    % Extract the data for this line
    Pts = C(:,pos+1:pos+Ni);
    % Plot the line
    [~,ind] = sort(dims);
    H = line(Pts(ind(1),:),Pts(ind(2),:),Pts(ind(3),:));
    set(H,'Color','w')
    % Re-position the pointer 
    pos = pos+Ni+1;
end    
