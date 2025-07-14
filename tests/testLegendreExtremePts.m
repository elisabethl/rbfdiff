function P=testLegendreExtremePts(l)
close all
%
% Test how the approximation behaves near -1,0,1
%
n = 500;
pts = logspace(-16,-1,n)';
pts0 = logspace(-64,-1,n)';
off = [-1 0 0 1];
mul = [1 -1 1 -1];
x{1} = -1+pts;
x{2} = -pts0;
x{3} = +pts0;
x{4} = 1-pts;
xstr = {'Distance from -1','Distance from $0^-$','Distance from $0^+$','Distance from 1'};
tstr = {'$|P_\mu^\nu|$', ...
        '$|\frac{\partial P_\mu^\nu}{\partial\theta}|$', ...
        '$|\frac{P_\mu^\nu}{\sin\theta}|$',...
        '$|\frac{\partial^2 P_\mu^\nu}{\partial\theta^2}|$', ...
        '$|A|$','$|B|$'};
for k = 1:size(x,2)
    [P{1},P{2},P{3},P{4},P{5},P{6}] = normalizedLegendre(l, x{k});
    for j=1:length(P)
        %
        % Make a figure for each type of function
        %
        if length(find(isnan(P{j})))>0
            j,k
            return
        end
        figure
        loglog(abs(x{k}-off(k)),abs(P{j}))
        adjustFigProp;
        xlabel(xstr{k},'interpreter','latex')
        title(tstr{j},'interpreter','latex')
    end
end

function adjustFigProp()
set(gca,'LineWidth',1.5,'FontSize',18)
set(gcf,'Color','w')
axis tight
grid

