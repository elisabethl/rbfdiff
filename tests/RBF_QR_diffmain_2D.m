clear

domain = 'square'; % The options are disc and square
N = 40; % The number of node points in the domain
Ne = 400; % The number of evaluation points
shift = 0; % To shift the center of the eval pts
scale = 0.8; % To scale the size of the evaluation area
epvec = logspace(-2,0,50); % To get error as a function of ep
epPlot = 0.1; % The value of ep for plotting the sol and error over the domain

fnum = 1;

% Add some functions to choose from and compute their derivatives symbolically
syms x y
if fnum==1
    f = 0.5*(x^0+y^0);
elseif fnum==2
    f = x + y;
elseif fnum==3
    f = x^2 + y^2;
elseif fnum==4
    f = sin(x-2*y);
elseif fnum==5
    f = 25/(1+25*(x^2+0.5*y^2));
elseif fnum==6
    f = sin(x*y+y^2)-cos(x^2);
end   
fun = matlabFunction(f,'vars',{'x','y'});
funx = matlabFunction(diff(f,x),'vars',{'x','y'});
funy = matlabFunction(diff(f,y),'vars',{'x','y'});
funxx = matlabFunction(diff(f,x,2),'vars',{'x','y'});
funxy = matlabFunction(diff(f,x,y),'vars',{'x','y'});
funyy = matlabFunction(diff(f,y,2),'vars',{'x','y'});


% Use Halton nodes in the domain, clustered towards the boundary
d=2;
if strcmp(domain,'square')
    xk = 2*(halton(N,d)-0.5);
    xk(:,1) = sin(pi/2*xk(:,1));
    xk(:,2) = sin(pi/2*xk(:,2));

    % Evaluation points
    ne = ceil(sqrt(Ne)); % Per dimension
    x = linspace(-1,1,ne);
    [xx,yy]=meshgrid(x);
    sz = size(xx);
    xe = [xx(:) yy(:)];
else % 'disc'
    Nin = ceil(1.2*4/pi*N);
    xk = 2*(halton(Nin,d)-0.5);
    [th,r]=cart2pol(xk(:,1),xk(:,2));
    pos = find(r<=1);
    pos = pos(1:N);
    xk = xk(pos,:); th = th(pos); r = r(pos);
    r = sin(pi/2*r);
    [xk(:,1),xk(:,2)] = pol2cart(th,r);

    % Evaluation points (Ne = nt*nr = 3*n*n)
    nr = ceil(sqrt(Ne/3));
    nt = 3*nr;
    th = 0:2*pi/nt:2*pi; % First and last is the same
    r = 0:1/(nr-1):1;
    [tt,rr] = meshgrid(th,r);
    sz = size(tt);
    [xe(:,1), xe(:,2)] = pol2cart(tt(:),rr(:));
end
xe = scale*xe+shift;

% Plot nodes
figure(1),clf
xx = reshape(xe(:,1)-shift,sz);
yy = reshape(xe(:,2)-shift,sz);
mesh(xx,yy,zeros(sz))
hold on
plot(xk(:,1),xk(:,2),'bo')

view(2)
axis equal
title('Node points: o, Evaluation points: grid')
%
% Apply to a function and compare with exact solution
%
uk = fun(xk(:,1),xk(:,2));
if length(uk)==1
    uk = uk*ones(size(xk(:,1)));
end    

epvec = [epvec  epPlot]
nep = length(epvec);

tic;
for j=1:length(epvec)
    ep = epvec(j);
    [A,Psi,Te] = RBF_QR_diffmat_2D(0,xk,ep,xe);
    [G,Te] = RBF_QR_diffmat_2D(1,Psi,Te);
    [L,Te] = RBF_QR_diffmat_2D(1.5,Psi,Te);
    [H,Te] = RBF_QR_diffmat_2D(2,Psi,Te);
    if (j<nep)
        maxerr(j) = max(abs(A*uk-fun(xe(:,1),xe(:,2))));
        errx(j) = max(abs(G{1}*uk-funx(xe(:,1),xe(:,2))));
        erry(j) = max(abs(G{2}*uk-funy(xe(:,1),xe(:,2))));
        errxx(j) = max(abs(H{1,1}*uk-funxx(xe(:,1),xe(:,2))));
        errxy(j) = max(abs(H{1,2}*uk-funxy(xe(:,1),xe(:,2))));
        erryy(j) = max(abs(H{2,2}*uk-funyy(xe(:,1),xe(:,2))));
        errL(j) = max(abs(L*uk-funxx(xe(:,1),xe(:,2))-funyy(xe(:,1),xe(:,2))));
    end     
end
epvec = epvec(1:end-1);
nep
time_nep=toc

figure(2)
loglog(epvec,maxerr,epvec,errx,epvec,erry,epvec,errxx,epvec,errxy,epvec,erryy,epvec,errL)

%
% Mesh the solution
%
figure(3)
ue = fun(xe(:,1),xe(:,2));
if length(ue)==1
    ue = ue*ones(size(xe(:,1)));
end    
mesh(xx,yy,reshape(ue,sz))
%
% Mesh the interpolation error, the x-derivative error, and the Laplacian error
%
figure(4)
y = A*uk-fun(xe(:,1),xe(:,2));
mesh(xx,yy,reshape(y,sz))

figure(5)
y = G{1}*uk-funx(xe(:,1),xe(:,2));             
mesh(xx,yy,reshape(y,sz))

figure(6)
y = L*uk-funxx(xe(:,1),xe(:,2))-funyy(xe(:,1),xe(:,2))
mesh(xx,yy,reshape(y,sz))
