close all
clear all
setPaths;

pars.dim = 2;                            % dim = 1,2 or 3
pars.display = 0;                        % Plot solution
pars.geom = 'ball';                      % ball or cube
pars.mode = 'unfitted';                  % fitted, unfitted or collocation
pars.bcMode = 'weak';                    % strong or weak imposition of boundary conditions (only relevant for fitted)
pars.scaling = 1;                        % Include scaling of the unfitted LS problem
pars.mvCentres = 1;                      % Option to have a Y point on top of all X points inside the domain
pars.q = 3;                              % Oversampling
                             % Number of center points (X) in each patch
pars.ep = 0.1;                           % Not relevant for 'r3' basis
pars.phi = 'rbfqr';                      % Choice of basis 'r3', 'mq', 'gs', 'iq', 'rbfqr'
pars.psi = 'wendland_c2';                % Weight function: wendland_c2 or bump
pars.pdeg = -1;                          % Polynomial extension, not relevant for 'rbfqr'
pars.del = 0.2;                          % Overlap between patches

order = 4;

P = [1 2 4 8 16 32 50];
P = P.^pars.dim;      
for i = 1:length(order)
    pars.N = nchoosek(order(i)+pars.dim,pars.dim);
    for j = 1:length(P)
        pars.P = P(j);                          % Number of patches
        pars.mode = "unfitted";
        [l2Error_RBFPU(i,j), h_RBFPU(i,j)] = RBFPUM_Poisson(pars);
        pars.mode = "collocation";
        [l2Error_RBFPU_C(i,j), h_RBFPU_C(i,j)] = RBFPUM_Poisson(pars);
    end
end

pars.rbfdeg = 4;
pars.extCoeff = 0.5;
pars.phi = 'r3';

N = [128 256 512 1024 2048 4096 8192 16384];
for i = 1:length(order)
    pars.pdeg = order(i);
    for j = 1:length(N)
        pars.N = N(j);     % Number of patches
        pars.mode = "unfitted";
        [l2Error_RBFFD(i,j), h_RBFFD(i,j)] = RBFFD_Poisson(pars);
        pars.mode = "collocation";
        [l2Error_RBFFD_C(i,j), h_RBFFD_C(i,j)] = RBFFD_Poisson(pars);
    end
end

figure(); 
for i = 1:length(order)
    loglog(1./h_RBFFD_C(i,:),l2Error_RBFFD_C(i,:),'r-s','LineWidth',1.5,'MarkerSize',12,'LineWidth',1.5)
    hold on
    A = [ones(length(l2Error_RBFFD_C),1), log(1./h_RBFFD_C')]; F = log(l2Error_RBFFD_C)'; slope = A\F;
    loglog(1./h_RBFFD_C,exp(slope(1)).*(1./h_RBFFD_C).^(slope(2)),'m--','LineWidth',2)

    loglog(1./h_RBFFD(i,:),l2Error_RBFFD(i,:),'r-s','LineWidth',1.5,'MarkerSize',12,'LineWidth',1.5)
    A = [ones(length(l2Error_RBFFD),1), log(1./h_RBFFD')]; F = log(l2Error_RBFFD)'; slope2 = A\F;
    loglog(1./h_RBFFD,exp(slope2(1)).*(1./h_RBFFD).^(slope2(2)),'m--','LineWidth',2)
    
    loglog(1./h_RBFPU_C(i,:),l2Error_RBFPU_C(i,:),'b-o','LineWidth',1.5,'MarkerSize',12,'LineWidth',1.5)
    A = [ones(length(l2Error_RBFPU_C),1), log(1./h_RBFPU_C')]; F = log(l2Error_RBFPU_C)'; slope3 = A\F;
    loglog(1./h_RBFPU_C,exp(slope3(1)).*(1./h_RBFPU_C).^(slope3(2)),'g--','LineWidth',2)

    loglog(1./h_RBFPU(i,:),l2Error_RBFPU(i,:),'b-o','LineWidth',1.5,'MarkerSize',12,'LineWidth',1.5)
    A = [ones(length(l2Error_RBFPU),1), log(1./h_RBFPU')]; F = log(l2Error_RBFPU)'; slope4 = A\F;
    loglog(1./h_RBFPU,exp(slope4(1)).*(1./h_RBFPU).^(slope4(2)),'g--','LineWidth',2)

    % loglog(1./hOutLS,l2ErrorLS,'b-o','LineWidth',1.5,'MarkerSize',12,'LineWidth',1.5)
    % hold on
    % A = [ones(length(l2ErrorLS),1), log(1./hOutLS')]; F = log(l2ErrorLS)'; slope2 = A\F;
    % loglog(1./hOutLS,exp(slope2(1)).*(1./hOutLS).^(slope2(2)),'g--','LineWidth',2)
    grid on;
    ax = gca;
    ax.FontSize = 18;
    xlabel("1/h","FontSize",28,"Interpreter","tex")
    ylabel("e_{rel}","FontSize",28,"Interpreter","tex",'Rotation',0)
    xlim([1, 100])
    ylim([5e-8 1])
    legend('Collocation RBF-FD',['slope = ',num2str(round(slope(2),2))],...
           'Unfitted LS-RBF-FD',['slope = ',num2str(round(slope2(2),2))],...
           'Collocation RBF-PUM',['slope = ',num2str(round(slope3(2),2))],...
           'Unfitted LS-RBF-PUM',['slope = ',num2str(round(slope4(2),2))],...
           "Interpreter","latex",'Location','south','NumColumns',2,'Orientation',...
           'horizontal','FontSize',18)

end
