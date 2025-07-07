function B = rbfmat_vpa_3D(ep,xk,xe,op)
digits(100);

ep = vpa(ep);

rbf = @(ep,rd2) exp(-ep^2*rd2);

[xi,xj] = meshgrid(xk(:,1)); xi = vpa(xi); xj = vpa(xj); 
[yi,yj] = meshgrid(xk(:,2)); yi = vpa(yi); yj = vpa(yj);
[zi,zj] = meshgrid(xk(:,3)); zi = vpa(zi); zj = vpa(zj);
% rd2 = (xi-xj).^2 + (yi-yj).^2 + (zi-zj).^2;
rd2 = vpa((xi-xj).^2 + (yi-yj).^2 + (zi-zj).^2);


% A = rbf(ep,rd2);
A = vpa(rbf(ep,rd2));
eigA = eig(A);
lambda = double(eigA);
disp(['The condition number of A is ' ...
    num2str(max(lambda)/min(lambda),'%.2e') '.']);

[xi,xj] = meshgrid(xk(:,1),xe(:,1)); xi = vpa(xi); xj = vpa(xj); 
[yi,yj] = meshgrid(xk(:,2),xe(:,2)); yi = vpa(yi); yj = vpa(yj);
[zi,zj] = meshgrid(xk(:,3),xe(:,3)); zi = vpa(zi); zj = vpa(zj);
% dx = xj-xi; dy = yj-yi; dz = zj-zi;
dx = vpa(xj-xi); dy = vpa(yj-yi); dz = vpa(zj-zi);
% rd2 = dx.^2 + dy.^2 + dz.^2;
rd2 = vpa(dx.^2 + dy.^2 + dz.^2);

if nargin<4
    Ae = rbf(ep,rd2);
else
    switch(op)
        case '0'
            Ae = rbf(ep,rd2);
        case 'x'
%             Ae = -2*ep^2*dx.*rbf(ep,rd2);
            Ae = vpa(-2*ep^2*dx.*rbf(ep,rd2));
        case 'y'
            Ae = -2*ep^2*dy.*rbf(ep,rd2);
        case 'z'
            Ae = -2*ep^2*dz.*rbf(ep,rd2);
        case 'L'
            Ae = 2*ep^2*(2*ep^2*rd2-3).*rbf(ep,rd2);
        case 'Ls'
            % Spherical Laplacian
            Ae = (-ep^4*rd2.^2+(4*ep^4+2*ep^2)*rd2-4*ep^2).*rbf(ep,rd2);
        otherwise
            disp('op not recognized.');
            return
    end
end

% B = double(Ae/A);
B = double(vpa(Ae/A));