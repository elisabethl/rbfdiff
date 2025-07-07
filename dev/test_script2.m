% %% Test script
% f = @(r) exp(-2*r.^2);
% Lf = @(r) 8*(8*r.^2-3).*f(r);
% n = 50;
% xk = 2*rand(n,3)-1;
% ep = 0.001;
% op = 'L1';
% delta = 1;
% % 
% % r = sqrt(xk(:,1).^2+xk(:,2).^2+xk(:,3).^2);
% % xk = xk./(delta*max(r));
% xk_c = xk;
% % 
% % rbf = @(ep,r2) exp(-ep^2*r2);
% % [xi,xj] = meshgrid(xk_c(:,1));
% % [yi,yj] = meshgrid(xk_c(:,2));
% % [zi,zj] = meshgrid(xk_c(:,3));
% % r2 = (xi-xj).^2+(yi-yj).^2+(zi-zj).^2;
% % 
% [la,th,r] = cart2sph(xk_c(:,1),xk_c(:,2),xk_c(:,3));
% % th = pi/2-th;
% xk_p = [r th la];
% % 
% % [Psi,Rt] = InitPsi_3D(ep,xk_p);
% 
% % ne = 5;
% % xe = 2*rand(ne,3)-1;
% xe = [0 0 0];
% % xe = xk;
% % r = sqrt(xe(:,1).^2+xe(:,2).^2+xe(:,3).^2);
% % if max(r)~=0
% %     xe = xe./(delta*max(r)); 
% % end
% % xe_c = xe;
% % [la,th,r] = cart2sph(xe_c(:,1),xe_c(:,2),xe_c(:,3));
% % th = pi/2-th;
% % xe_p = [r th la];
% 
% % A = RBF_QR_mat_3D(Psi,'0',xk_p);
% % Ae = RBF_QR_mat_3D(Psi,op,xe_p);
% % B = Ae/A; % RBF-QR Matrix for op xk->xe
% % Bx = RBF_QR_diffmat_3D(op,xe,xk,ep);
% 
% B = RBF_QR_diffmat_3D(op,xe,xk,ep);
% 
% % A_d = rbf(ep,r2);
% % [xi,xj] = meshgrid(xk_c(:,1),xe_c(:,1));
% % [yi,yj] = meshgrid(xk_c(:,2),xe_c(:,2));
% % [zi,zj] = meshgrid(xk_c(:,3),xe_c(:,3));
% % r2 = (xi-xj).^2+(yi-yj).^2+(zi-zj).^2;
% % Ae_d = rbf(ep,r2);
% % B_d = Ae_d/A_d; % RBF-Direct Matrix for interpolation xk->xe
% % B_v = rbfmat_vpa_3D(ep,xk_c,xe_c,op); %RBF-VPA matrix for op xk->xe
% 
% % B_v = rbfmat_vpa_3D(ep,xk,xe,op);
% % B_v = rbfmat_vpa_3D(ep,xk,xe,'L');
% 
% eptilde = ep*(1+1000*randn(1)*eps);
% B_v = RBF_QR_diffmat_3D(op,xe,xk,eptilde);
% 
% 
% fe = B*f(xk_p(:,1));
% % fe_d = B_d*f(xk_p(:,1));
% fe_v = B_v*f(xk_p(:,1));
% % 

% % % % %%% NEW TEST, diffmatrices on the sphere.
n = 400;
ep = 0.001;
delta = 1e-10; % Relative perturbation of epsilon
op = {'x','y','z'};
% op = 'L';
% op = 'x';
fd = 7;
nodes = load(['me' num2str(n) '.dat']);
x = nodes(:,1); y = nodes(:,2); z = nodes(:,3);
Dx = zeros(n);
Dy = zeros(n);
Dz = zeros(n);
% % Dx_v = zeros(n);
% % Dy_v = zeros(n);
% % Dz_v = zeros(n);
% D = zeros(n);
D_v = zeros(n);
for i = 1:n
% for i = 245;
    xe = [x(i) y(i) z(i)];
    r2 = (x-xe(1)).^2+(y-xe(2)).^2+(z-xe(3)).^2;
    [~,idx] = sort(r2,'ascend');
    idx_i = idx(1:fd);
    xk = [x(idx_i) y(idx_i) z(idx_i)];
    
    cc = sum(xk,1)/size(xk,1);
    xk_c = [xk(:,1)-cc(1) xk(:,2)-cc(2) xk(:,3)-cc(3)];
    r = sqrt(sum(xk_c.^2,2));
    rr = max(r);
    xk_c = xk_c/rr;
    [la,th,r] = cart2sph(xk_c(:,1),xk_c(:,2),xk_c(:,3));
    xk_p = [r pi/2-th la];
    xe_c = [xe(:,1)-cc(1) xe(:,2)-cc(2) xe(:,3)-cc(3)]/rr;
    [la,th,r] = cart2sph(xe_c(:,1),xe_c(:,2),xe_c(:,3));
    xe_p = [r pi/2-th la];
%     Psi = InitPsi_3D(ep*rr,xk_p);
%     A = RBF_QR_mat_3D(Psi,'0',xk_p);
%     Ae = RBF_QR_mat_3D(Psi,op,xe_p);
%     W1 = Ae/A;
%     W1 = W1/rr;
    
    [W,Psi] = RBF_QR_diffmat_3D(op,xe,xk,ep);
%     D(i,idx_i) = W;
    Dx(i,idx_i) = W{1};
    Dy(i,idx_i) = W{2};
    Dz(i,idx_i) = W{3};
    
%     eptilde = ep*(1+delta*randn(1));
%     W_v = RBF_QR_diffmat_3D('L',xe,xk,ep);
%     W_v = RBF_QR_diffmat_3D(op,xe,xk,eptilde);
%     W_v = rbfmat_vpa_3D(ep,xk,xe,op);
%     W_v = rbfmat_vpa_3D(ep,xk,xe,'L');
    W_v = rbfmat_vpa_3D(ep,xk,xe,'Ls');
%     W_v = rbfmat_vpa_3D(eptilde,xk,xe,'L');
%     W_v = rbfmat_vpa_3D(ep*rr,xk_c,xe_c,op)/rr^2;
    D_v(i,idx_i) = W_v;
%     rel_diff = abs(W-W_v)./abs(W_v);
%     W_v_ = rbfmat_vpa_3D(eptilde,xk,xe,op{1});
%     W_v = rbfmat_vpa_3D(ep,xk,xe,op{1});
%     W_v = rbfmat_vpa_3D(ep*rr,xk_c,xe_c,op{1})/rr;
%     Dx_v(i,idx_i) = W_v;
%     W_v = rbfmat_vpa_3D(ep*rr,xk_c,xe_c,op{2})/rr;
%     Dy_v(i,idx_i) = W_v;
%     W_v = rbfmat_vpa_3D(ep*rr,xk_c,xe_c,op{3})/rr;
%     Dz_v(i,idx_i) = W_v;
    rel_diff = abs(W{1}-W_v)./abs(W_v);
    disp(num2str([i norm(rel_diff)]));
%     disp(num2str([max(abs(W_v_-W_v)) max(abs(W_v)) max(abs(W_v_))]));
%     disp(num2str([max(abs(W{1}-W_v)) max(abs(W_v)) max(abs(W{1}))]));
%     pause();
%     rel_diff2 = norm(abs(W_v-W_v_)./abs(W_v_))
%     if norm(rel_diff)>1e-8
% %         disp(['Stencil ' num2str(i) ': ' num2str(norm(rel_diff)) '.']);
% % %         [xx,yy,zz] = sph2cart(Psi.xk(:,3),pi/2-Psi.xk(:,2),Psi.xk(:,1));
% % %         scatter3(xx,yy,zz,50,rel_diff,'filled');
% % %         pause(0.1);
% %         disp(['Stencil ' num2str(i) ' at ' num2str(xe) ':']);
%         disp(['Weights (QR): ' num2str(W)]);
% %         disp(['Weights (.1): ' num2str(W1)]);
%         disp(['Weights (VP): ' num2str(W_v)]);
%     end
end

% % [DPx_,DPy_,DPz_,L] = rbfmatrix_fd_hyper([x y z],ep,fd,1,3);
% % DPx_ = full(DPx_);
% % [DPx,DPy,DPz] = dmcart2proj(x,y,z,Dx,Dy,Dz);
% % [DPx_v,DPy_v,DPz_v] = dmcart2proj(x,y,z,Dx_v,Dy_v,Dz_v);
