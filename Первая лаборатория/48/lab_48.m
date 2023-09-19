load('data_48');
points = [1,2,3,4,6]; % Номера параметров, участвующих в расчёте средних значений

[w, Ra, mu_a] = deal(zeros(size(B)));
fit_H_B = cell(size(B));
for i = 1:numel(B)
    w(i) = trapz(B{i},H{i}); % Дж/м^3
    % Аппроксимация B(H)
    fit_H_B{i} = fit(H{i}(1:50),B{i}(1:50),{'x','x^2'});
    [mu_a(i), Ra(i)] = matsplit(coeffvalues(fit_H_B{i})./mu_0);
end
% Средние значения параметров
[mean_mu_a, d_mean_mu_a] = samplemean(mu_a(points));
[mean_Ra, d_mean_Ra] = samplemean(Ra(points)); % м/А
[mean_B_s, d_mean_B_s] = samplemean(B_s(points)); % Тл
[mean_B_r, d_mean_B_r] = samplemean(B_r(points)); % Тл
[mean_H_c, d_mean_H_c] = samplemean(H_c(points)); % А/м
[mean_w, d_mean_w] = samplemean(w(points)); % Дж/м^3
%% Вывод результатов
printresult(mu_a,4,'s','μ_a');
printresult(mean_mu_a,d_mean_mu_a,'pm','mean μ_a');
printresult(Ra,4,'s','Ra','м/А');
printresult(mean_Ra,d_mean_Ra,'pm','mean Ra','м/А');
printresult(B_s,4,'s','B_s','Тл');
printresult(mean_B_s,d_mean_B_s,'pm','mean B_s','Тл');
printresult(B_r,3,'s','B_r','Тл');
printresult(mean_B_r,d_mean_B_r,'pm','mean B_r','Тл');
printresult(H_c,4,'s','H_c','А/м');
printresult(mean_H_c,d_mean_H_c,'pm','mean H_c','А/м');
printresult(w,4,'s','w','Дж/м^3');
printresult(mean_w,d_mean_w,'pm','mean w','Дж/м^3');
%% B(H)
i = 1;

figure('Name','B(H)','Units','centimeters');
plot(H{i},B{i},'Color','k');
set(gca,'FontName','Calibri','FontSize',13,'LooseInset',get(gca,'TightInset'));
xlabel('\itH\rm, А/м','FontName','Liberation Serif');
ylabel('\itB\rm, Тл','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',10);