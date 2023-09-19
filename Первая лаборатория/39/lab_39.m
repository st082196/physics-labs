load('data_39');
spectrum = imread('img\visible_spectrum_380-750.png');

% Спектры совместного пропускания
T.red_yellow.theor = T.red.*T.yellow.*1e-2;
T.green_turquoise.theor = T.green.*T.turquoise.*1e-2;

% Аппроксимация C(E)
fit_E_C = cell(numel(fit_points),1);
C_unknown = zeros(numel(fit_points),1);
for i = 1:numel(fit_points)
    fit_E_C{i} = fit(E(fit_points(i),:).',C.',{'x'});
    C_unknown(i) = fit_E_C{i}(E_unknown(fit_points(i)));
end
% Концентрация неизвестного раствора, %
mean_C_unknown = mean(C_unknown);
d_mean_C_unknown = sqrt(var(C_unknown)./numel(C_unknown));

% ∆_C, %
fit_E_C_shifted = cell(numel(fit_points),numel(C));
C_unknown_shifted = zeros(numel(fit_points),numel(C));
for j = 1:numel(C)
    C_shifted = C;
    C_shifted(j) = C_shifted(j) + d_C(j);
    for i = 1:numel(fit_points)
        fit_E_C_shifted{i,j} = fit(E(fit_points(i),:).',C_shifted.',{'x'});
        C_unknown_shifted(i,j) = fit_E_C_shifted{i,j}(E_unknown(fit_points(i)));
    end
end
mean_C_unknown_shifted = mean(C_unknown_shifted);
d_C_unknown = norm(mean_C_unknown_shifted - mean_C_unknown);
%% Вывод результатов
printresult(mean_C_unknown,d_mean_C_unknown,'pm','C','%',-2);
printresult(d_C_unknown,2,'s','∆_C','%',-2);
%% I_ref(λ)
figure('Name','I_ref(λ)','Units','centimeters');
hold on;

image('CData',spectrum,'XData',[380 750],'YData',[0 100]);
plot(wl,I.water,'Color','w');

set(gca,'FontName','Calibri','FontSize',14,'Color','k','Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itI\rm_{ref}, %','FontName','Liberation Serif');
xlim([407 996]);
ylim([0 100]);
pbaspect([2 1 1]);
figresize('width',17);
%% T(λ), отдельный фильтр
color = 'red';
[titles.red, titles.yellow, titles.green, titles.turquoise] = ...
    matsplit(["Красный фильтр", "Жёлтый фильтр", "Зелёный фильтр", "Бирюзовый фильтр"]);

figure('Name','T(λ)','Units','centimeters');
hold on;

image('CData',spectrum,'XData',[380 750],'YData',[0 100]);
plot(wl,T.(color),'Color','w');

set(gca,'FontName','Calibri','FontSize',14,'Color','k','Box','on','LooseInset',get(gca,'TightInset'));
title(titles.(color));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itT\rm, %','FontName','Liberation Serif');
xlim([407 996]);
ylim([0 100]);
pbaspect([2 1 1]);
figresize('width',17);
%% T(λ), комбинация фильтров
color = 'red_yellow';
[titles.red_yellow, titles.green_turquoise] = ...
    matsplit(["Красный + жёлтый фильтр", "Зелёный + бирюзовый фильтр"]);

figure('Name','T(λ)','Units','centimeters');
hold on;

image('CData',spectrum,'XData',[380 750],'YData',[0 100]);
plot(wl,T.(color).exp,'Color','w');
plot(wl,T.(color).theor,'Color','b');

set(gca,'FontName','Calibri','FontSize',14,'Color','k','Box','on','LooseInset',get(gca,'TightInset'));
title(titles.(color));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itT\rm, %','FontName','Liberation Serif');
legend('измеренный спектр','рассчитанный спектр','Location','southeast','Box','off','TextColor','w');
xlim([407 996]);
ylim([0 100]);
pbaspect([2 1 1]);
figresize('width',17);
%% E(λ)
figure('Name','E(λ)','Units','centimeters');
hold on;

image('CData',spectrum,'XData',[380 750],'YData',[0 1]);
color = ["b","g","m"];
for j = 1:numel(color)
    plot(wl,E(:,j),'Color',color(j));
end
plot(wl,E_unknown,'Color','w');

set(gca,'FontName','Calibri','FontSize',14,'Color','k','Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itE','FontName','Liberation Serif');
legend('раствор №1','раствор №2','раствор №3','неизвестный раствор','Location','northeast','Box','off','TextColor','w');
xlim([407 996]);
ylim([0 1]);
pbaspect([2 1 1]);
figresize('width',17);
%% C(E)
i = 736;

figure('Name','C(E)','Units','centimeters');
hold on;

x_fit = [0 1];
plot(x_fit,fit_E_C{fit_points == i}(x_fit),'Color','g');
plot(E(i,:),C,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(E_unknown(i,:),C_unknown(fit_points == i),'o','Color','b','MarkerFaceColor','w','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title("\itλ\rm\bf = " + wl(i) + " нм",'FontName','Liberation Serif');
xlabel('\itE','FontName','Liberation Serif');
ylabel('\itC\rm, %','FontName','Liberation Serif');
legend('аппроксимирующая прямая','известные растворы','неизвестный раствор','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% C(λ)
figure('Name','C(λ)','Units','centimeters');
hold on;

image('CData',spectrum,'XData',[380 750],'YData',[0.014 0.016]);
plot(wl(fit_points),C_unknown,'.','Color','w');
yline(mean_C_unknown,'Color','b');

set(gca,'FontName','Calibri','FontSize',14,'Color','k','Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itC\rm, %','FontName','Liberation Serif');
legend('спектр значений концентрации','среднее значение','Location','northeast','TextColor','w','Color','none');
xlim(wl(fit_points([1 end])));
ylim([0.014 0.016]);
pbaspect([2 1 1]);
figresize('width',17);