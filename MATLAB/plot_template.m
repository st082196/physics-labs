figure('Name','y(x)','Units','centimeters');
hold on;

X = [0 1];
marker = {'o','s','^','v','d','p','h'};
color = {'b','r','g','m','c','#A2142F','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(X,fitobject{j}(X),'Color',color{j});
    subset(j) = plot(x{j},y{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, А','FontName','Liberation Serif');
ylabel('\ity\rm, В','FontName','Liberation Serif');
legend(subset,'№1','№2','№3','№4','№5','Location','southeast');
% xlim([0 1]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);