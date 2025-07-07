clear
close all

domain = 'sphere'; % The options are sphere and cube
N = 40; % The number of node points in the domain
Ne = 1000; % The number of evaluation points
shift = 0; % To shift the center of the eval pts
scale = 0.8; % To scale the size of the evaluation area
epvec = logspace(-2,0,25); % To get error as a function of ep
epPlot = 0.1;
rPlot = [0.5 0.8 1]; % Radii to plot for or z-values

fnum = 1;

% Add some functions to choose from and compute their derivatives symbolically
syms x y z
if fnum==1
    f = 0.5*(x^0+y^0+z^0);
elseif fnum==2
    f = x + y + z;
elseif fnum==3
    f = x^2 + y^2 - z^2;
elseif fnum==4
    f = sin(x-2*y+0.5*z);
elseif fnum==5
    f = 25/(1+25*(x^2+0.5*y^2+0.1*z^2));
elseif fnum==6
    f = sin(x*y+y^2)-cos(x^2-z^2);
end   
fun = matlabFunction(f,'vars',{'x','y','z'});
funx = matlabFunction(diff(f,x),'vars',{'x','y','z'});
funy = matlabFunction(diff(f,y),'vars',{'x','y','z'});
funz = matlabFunction(diff(f,z),'vars',{'x','y','z'});
funxx = matlabFunction(diff(f,x,2),'vars',{'x','y','z'});
funxy = matlabFunction(diff(f,x,y),'vars',{'x','y','z'});
funxz = matlabFunction(diff(f,x,z),'vars',{'x','y','z'});
funyy = matlabFunction(diff(f,y,2),'vars',{'x','y','z'});
funyz = matlabFunction(diff(f,y,z),'vars',{'x','y','z'});
funzz = matlabFunction(diff(f,z,2),'vars',{'x','y','z'});

% Use Halton nodes in the domain, clustered towards the boundary
d=3;
if strcmp(domain,'cube')
    xk = 2*(halton(N,d)-0.5);
    xk(:,1) = sin(pi/2*xk(:,1));
    xk(:,2) = sin(pi/2*xk(:,2));
    xk(:,3) = sin(pi/2*xk(:,3));

    % Evaluation points
    ne = ceil(Ne^(1/d)); % Per dimension
    x = linspace(-1,1,ne);
    [xx,yy,zz]=ndgrid(x);
    sz = size(xx);
    xe = [xx(:) yy(:) zz(:)];
else % 'sphere'
    Nin = ceil(1.2*8/(4*pi/3)*N);
    xk = 2*(halton(Nin,d)-0.5);
    [th,phi,r]=cart2sph(xk(:,1),xk(:,2),xk(:,3));
    pos = find(r<=1);
    pos = pos(1:N);
    xk = xk(pos,:); th = th(pos); phi = phi(pos); r = r(pos);
    r = sin(pi/2*r);
    [xk(:,1),xk(:,2),xk(:,3)] = sph2cart(th,phi,r);

    % Evaluation points (Ne = nsph*nr = 2*n^2*n)
    nr = ceil((Ne/2)^(1/3));
    nsph = floor(sqrt(2)*nr);
    [xx,yy,zz] = sphere(nsph);
    sz = size(xx);
    r = 0:1/(nr-1):1;
    xe = zeros(0,d);
    for k=1:length(r)
        xe = [xe; r(k)*[xx(:) yy(:) zz(:)]];
    end
    sz = [sz length(r)];
end
xe = scale*xe+shift;

% Plot nodes
figure(1),clf
plot3(xk(:,1),xk(:,2),xk(:,3),'bo')
axis equal
title('Node points')
%
% Apply to a function and compare with exact solution
%
uk = fun(xk(:,1),xk(:,2));
if length(uk)==1
    uk = uk*ones(size(xk(:,1)));
end    

epvec = [epvec  epPlot];
nep = length(epvec);

tic;
for j=1:length(epvec)
    ep = epvec(j);
    [A,Psi,Te] = RBF_QR_diffmat_3D(0,xk,ep,xe);
    [G,Te] = RBF_QR_diffmat_3D(1,Psi,Te);
    [L,Te] = RBF_QR_diffmat_3D(1.5,Psi,Te);
    [H,Te] = RBF_QR_diffmat_3D(2,Psi,Te);
    if (j<nep)
        maxerr(j) = max(abs(A*uk-fun(xe(:,1),xe(:,2),xe(:,3))));
        errx(j) = max(abs(G{1}*uk-funx(xe(:,1),xe(:,2),xe(:,3))));
        erry(j) = max(abs(G{2}*uk-funy(xe(:,1),xe(:,2),xe(:,3))));
        errz(j) = max(abs(G{2}*uk-funz(xe(:,1),xe(:,2),xe(:,3))));
        errxx(j) = max(abs(H{1,1}*uk-funxx(xe(:,1),xe(:,2),xe(:,3))));
        errxy(j) = max(abs(H{1,2}*uk-funxy(xe(:,1),xe(:,2),xe(:,3))));
        errxz(j) = max(abs(H{1,3}*uk-funxz(xe(:,1),xe(:,2),xe(:,3))));
        erryy(j) = max(abs(H{2,2}*uk-funyy(xe(:,1),xe(:,2),xe(:,3))));
        erryz(j) = max(abs(H{2,3}*uk-funyz(xe(:,1),xe(:,2),xe(:,3))));
        errzz(j) = max(abs(H{3,3}*uk-funzz(xe(:,1),xe(:,2),xe(:,3))));
        errL(j) = max(abs(L*uk-funxx(xe(:,1),xe(:,2))-funyy(xe(:,1),xe(:,2))));
    end     
end
epvec = epvec(1:end-1);
nep
time_nep=toc

figure(2)
loglog(epvec,maxerr,epvec,errx,epvec,erry,epvec,errz, ...
       epvec,errxx,epvec,errxy,epvec,errxz,epvec,erryy, ...
       epvec,erryz,epvec,errzz,epvec,errL)

%
% Mesh the interpolation error, the y-derivative error, and the Laplacian error
%
for k=1:length(rPlot)
    if (strcmp(domain,'sphere'))
        vec = r;
    else
        vec = x;
    end
    [~,pos] = min(abs(vec-rPlot(k)));
    trueval = vec(pos);

    
    figure
    y = A*uk-fun(xe(:,1),xe(:,2),xe(:,3));
    plotErr(domain,xx,yy,zz,y,sz,pos)
    title(['Interpolation error at xyz or r = ' num2str(trueval) ...
           ' for ep = ' num2str(epPlot)])

    figure
    y = G{2}*uk-funy(xe(:,1),xe(:,2),xe(:,3));             
    plotErr(domain,xx,yy,zz,y,sz,pos)
    title(['y-derivative error at xyz or r = ' num2str(trueval)])

    figure
    y = L*uk-funxx(xe(:,1),xe(:,2),xe(:,3))-funyy(xe(:,1),xe(:,2),xe(:,3)) - ...
        funzz(xe(:,1),xe(:,2),xe(:,3));
    plotErr(domain,xx,yy,zz,y,sz,pos)
    title(['Laplacian error at xyz or r = ' num2str(trueval)])
end
    
function str = plotErr(domain,xx,yy,zz,u,sz,pos,trueval)
    uu = log10(abs(reshape(u,sz)));
    if (strcmp(domain,'sphere'))
        H=surf(xx,yy,zz,uu(:,:,pos));
        set(H,'EdgeColor','none','FaceColor','interp')
    else
        H = mesh(xx(:,:,pos),yy(:,:,pos),zz(:,:,pos),uu(:,:,pos));
        set(H,'EdgeColor','none','FaceColor','interp')
        hold on
        x = squeeze(xx(:,pos,:)); y = squeeze(yy(:,pos,:));
        z = squeeze(zz(:,pos,:)); u = squeeze(uu(:,pos,:));
        H = mesh(x,y,z,u);
        set(H,'EdgeColor','none','FaceColor','interp')
        x = squeeze(xx(pos,:,:)); y = squeeze(yy(pos,:,:));
        z = squeeze(zz(pos,:,:)); u = squeeze(uu(pos,:,:));
        H = mesh(x,y,z,u);
        set(H,'EdgeColor','none','FaceColor','interp')
    end
    axis equal
    colorbar
end

