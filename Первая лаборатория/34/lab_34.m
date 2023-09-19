load('data_34');
%% 2.
unpackstruct(part2);
%% I_а(U_а)
figure('Name','I_а(U_а)','Units','centimeters');
hold on;

marker = {'o','s','^','v'};
color = {'b','r','g','m'};
for j = numel(color):-1:1
    plot(U_a,I_a(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itU\\rm_{нак} = %.1f В',U_h(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_а, В','FontName','Liberation Serif');
ylabel('\itI\rm_а, мА','FontName','Liberation Serif');
legend('Location','northwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% I_а(U_а^(3/2))
figure('Name','I_а(U_а^(3/2))','Units','centimeters');
hold on;

marker = {'o','s','^','v'};
color = {'b','r','g','m'};
for j = numel(color):-1:1
    plot(U_a.^(3./2),I_a(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itU\\rm_{нак} = %.1f В',U_h(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_а^{3/2}, В^{3/2}','FontName','Liberation Serif');
ylabel('\itI\rm_а, мА','FontName','Liberation Serif');
legend('Location','northwest','FontName','Liberation Serif');
xlim([0 1000]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.1.
unpackstruct(part3(1));
%% I_а(U_а)
figure('Name','I_а(U_а)','Units','centimeters');
hold on;

marker = {'o','s','^','v','d'};
color = {'b','r','g','m','k'};
for j = numel(color):-1:1
    plot(U_a,I_a(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itU\\rm_с = %.0f В',U_g(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_а, В','FontName','Liberation Serif');
ylabel('\itI\rm_а, мА','FontName','Liberation Serif');
legend('Location','northwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.2.
unpackstruct(part3(2));
%% I_а(U_с)
figure('Name','I_а(U_с)','Units','centimeters');
hold on;

marker = {'o','s','^'};
color = {'b','r','g'};
for j = numel(color):-1:1
    plot(U_g,I_a(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itU\\rm_а = %.0f В',U_a(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_с, В','FontName','Liberation Serif');
ylabel('\itI\rm_а, мА','FontName','Liberation Serif');
legend('Location','northwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);