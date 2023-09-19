load('data_12'); %#ok<*SAGROW>

[a, d_a, b, d_b] = deal(zeros(size(fit1)));
for j = 1:numel(fit1)
    fit1(j).points = r{j}.^-2 >= fit1(j).range(1) & r{j}.^-2 <= fit1(j).range(2);
    fit1(j).x = r{j}(fit1(j).points).^-2;
    fit1(j).y = F{j}(fit1(j).points);
    fit1(j).f = fit(fit1(j).x,fit1(j).y,{'x','1'});
    [a(j), b(j)] = matsplit(coeffvalues(fit1(j).f)); % мН⋅см^2, мН
    [d_a(j), d_b(j)] = matsplit(coeffvariances(fit1(j).f));
end
% U_C, В
[U_C, d_U_C] = samplemean(sample_U_C,d_U_C_syst);
% q, нКл
q = C.*U_C;
d_q = sqrt((U_C.*d_C).^2 + (C.*d_U_C).^2);
% Аппроксимация q(U)
fit2 = fit(U,q,{'x'});
% C_сф, пФ
C_sph.exp = coeffvalues(fit2);
d_C_sph.exp = coeffvariances(fit2);
C_sph.theor = 4.*pi.*e_0.theor.*R.*1e10;
d_C_sph.theor = 4.*pi.*e_0.theor.*d_R.*1e10;
% Аппроксимация a(q)
fit3.x = q(fit3.x_points);
fit3.y = a(fit3.y_points).';
fit3.f = fit(fit3.x,fit3.y,{'x'});
c = coeffvalues(fit3.f).*1e2; % Н⋅м^2⋅Кл^−1
d_c = coeffvariances(fit3.f).*1e2;
% Q, нКл
Q = q(1);
d_Q = d_q(1);
% ε_0, Ф/м
e_0.exp = Q./(4.*pi.*c).*1e-9;
d_e_0.exp = abs(e_0.exp).*sqrt((d_Q./Q).^2 + (d_c./c).^2);
%% Вывод результатов
printresult(a,d_a,'pm','a','мН⋅см^2');
printresult(b,d_b,'pm','b','мН');
printresult(U_C,d_U_C,'pm','U_C','В');
printresult(q,d_q,'pm','q','нКл');
printresult(C_sph.exp,d_C_sph.exp,'pm','C_сф, эксп','пФ');
printresult(C_sph.theor,d_C_sph.theor,'pm','C_сф, теор','пФ');
printresult(c,d_c,'pm','c','Н⋅м^2⋅Кл^−1');
printresult(e_0.exp,d_e_0.exp,'pm','ε_0, эксп','Ф/м',-12);
%% F(r^−2)
figure('Name','F(r^−2)','Units','centimeters');
hold on;

marker = {'o','s','^','v','d','p','h'};
color = {'b','r','g','m','c','#A2142F','k'};
for j = 1:numel(color)
    plot(r{j}.^-2,F{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itr\rm^{−2}, см^{−2}','FontName','Liberation Serif');
ylabel('\itF\rm, мН','FontName','Liberation Serif');
legend( ...
    '\itU\rm = +25.0 кВ', ...
    '\itU\rm = +15.0 кВ', ...
    '\itU\rm = +5.0 кВ', ...
    '\itU\rm = 0.0 кВ', ...
    '\itU\rm = −5.0 кВ', ...
    '\itU\rm = −15.0 кВ', ...
    '\itU\rm = −25.0 кВ', ...
    'Location','southwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% F(r^−2), линейный участок
figure('Name','F(r^−2)','Units','centimeters');
hold on;

marker = {'o','s','^','v','d','p','h'};
color = {'b','r','g','m','c','#A2142F','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(fit1(j).range,fit1(j).f(fit1(j).range),'Color',color{j});
    subset(j) = plot(fit1(j).x,fit1(j).y,marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('Линейный участок');
xlabel('\itr\rm^{−2}, см^{−2}','FontName','Liberation Serif');
ylabel('\itF\rm, мН','FontName','Liberation Serif');
legend(subset, ...
    '\itU\rm = +25.0 кВ', ...
    '\itU\rm = +15.0 кВ', ...
    '\itU\rm = +5.0 кВ', ...
    '\itU\rm = 0.0 кВ', ...
    '\itU\rm = −5.0 кВ', ...
    '\itU\rm = −15.0 кВ', ...
    '\itU\rm = −25.0 кВ', ...
    'Location','southwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% q(U)
figure('Name','q(U)','Units','centimeters');
hold on;

x_fit = [-25 25];
plot(x_fit,fit2(x_fit),'Color','g');
plot(U,q,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm, кВ','FontName','Liberation Serif');
ylabel('\itq\rm, нКл','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% a(q)
figure('Name','a(q)','Units','centimeters');
hold on;

x_fit = [-50 50];
plot(x_fit,fit3.f(x_fit),'Color','g');
plot(fit3.x,fit3.y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(q(end),a(end),'x','Color','r','MarkerFaceColor','r','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itq\rm, нКл','FontName','Liberation Serif');
ylabel('\ita\rm, мН·см^2','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);