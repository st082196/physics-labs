load('data_25');
%% 1.
unpackstruct(part1);

% Аппроксимация E(U)
[x, y, fit1] = deal(cell(size(E)));
[a, d_a, b, d_b] = deal(zeros(size(fit1)));
for j = 1:numel(fit1)
    x{j} = U{j}(U{j} < 10);
    y{j} = E{j}(U{j} < 10);
    fit1{j} = fit(x{j},y{j},{'x','1'});
    [a(j), b(j)] = matsplit(coeffvalues(fit1{j})); % м^−1, кВ/м
    [d_a(j), d_b(j)] = matsplit(coeffvariances(fit1{j}));
end
% d, мм
d = 1./a.*1e3;
d_d = d_a./a.^2.*1e3;
%% Вывод результатов
printresult(a,d_a,'pm','a','м^−1');
printresult(b,d_b,'pm','b','кВ/м');
printresult(d,d_d,'pm','d','мм');
%% E(U)
figure('Name','E(U)','Units','centimeters');
hold on;

x_fit = [0 10];
marker = {'o','s','^','v'};
color = {'b','r','g','m'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(x_fit,fit1{j}(x_fit),'Color',color{j});
    subset(j) = plot(x{j},y{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm, кВ','FontName','Liberation Serif');
ylabel('\itE\rm, кВ/м','FontName','Liberation Serif');
legend(subset, ...
    '\itd\rm = 10 мм','\itd\rm = 20 мм','\itd\rm = 100 мм','\itd\rm = 300 мм', ...
    'Location','northwest','FontName','Liberation Serif');
ylim([0 800]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% Аппроксимация E(d^−1)
fit2.points = d.^-1 >= fit2.range(1) & d.^-1 <= fit2.range(2);
fit2.x = d(fit2.points).^-1;
fit2.y = E(fit2.points);
fit2.f = fit(fit2.x,fit2.y,{'x'});
fit2.U = coeffvalues(fit2.f).*1e-3; % кВ
fit2.d_U = coeffvariances(fit2.f).*1e-3;
% Аппроксимация E(d)
fit3.points = d >= fit3.range(1) & d <= fit3.range(2);
fit3.x = d(fit3.points);
fit3.y = E(fit3.points);
fit3.f = fit(fit3.x,fit3.y,'U/abs(x-d_0)*1e3','StartPoint',[1 0]);
[fit3.U, fit3.d_0] = matsplit(coeffvalues(fit3.f)); % кВ, мм
[fit3.d_U, fit3.d_d_0] = matsplit(coeffvariances(fit3.f));
%% Вывод результатов
fprintf('\nЛинейная аппроксимация\n');
printresult(fit2.U,fit2.d_U,'pm','U','кВ');
fprintf('\nНелинейная аппроксимация\n');
printresult(fit3.U,fit3.d_U,'pm','U','кВ');
printresult(fit3.d_0,fit3.d_d_0,'pm','d_0','мм');
%% E(d)
figure('Name','E(d)','Units','centimeters');
plot(d,E,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itd\rm, мм','FontName','Liberation Serif');
ylabel('\itE\rm, кВ/м','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% E(d^−1)
figure('Name','E(d^−1)','Units','centimeters');
hold on;

plot(fit2.range,fit2.f(fit2.range),'Color','g');
plot(fit2.x,fit2.y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itd\rm^{−1}, мм^{−1}','FontName','Liberation Serif');
ylabel('\itE\rm, кВ/м','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% E(d), нелинейная аппроксимация
figure('Name','E(d)','Units','centimeters');
hold on;

x_fit = linspace(fit3.range(1),fit3.range(2),500);
plot(x_fit,fit3.f(x_fit),'Color','g');
plot(fit3.x,fit3.y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itd\rm, мм','FontName','Liberation Serif');
ylabel('\itE\rm, кВ/м','FontName','Liberation Serif');
ylim([0 600]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);