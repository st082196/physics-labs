%% Диффузионная длина
load('data_90');

% Аппроксимация ln(V/V_0) = f(x)
fit1.f = fit(x,log(V./V(1)),{'x'});
fit1.a = coeffvalues(fit1.f); % мм^−1
fit1.d_a = coeffvariances(fit1.f);
fit1.l_0 = -1./fit1.a; % мм
fit1.d_l_0 = abs(fit1.d_a./fit1.a.^2);
% Аппроксимация V(x)
[fit2.f, fit2.gof] = fit(x,V,'V0.*exp(-sqrt((x-x0).^2+y0.^2)./l0)', ...
    'StartPoint',[300,1,0,0],'Robust','LAR','MaxFunEvals',1e4,'MaxIter',1e4);
[fit2.V_0, fit2.l_0, fit2.x_0, fit2.y_0] = matsplit(coeffvalues(fit2.f)); % мВ, мм, мм, мм
[fit2.d_V_0, fit2.d_l_0, fit2.d_x_0, fit2.d_y_0] = matsplit(coeffvariances(fit2.f));
%% Вывод результатов
printresult(x,2,'d','x','мм');
printresult(V,0,'d','V','мВ');
printresult(fit1.a,fit1.d_a,'pm','a','мм^−1');
printresult(fit1.l_0,fit1.d_l_0,'pm','l_0','мм');
printresult(fit2.V_0,fit2.d_V_0,'pm','V_0','мВ',3);
printresult(fit2.l_0,fit2.d_l_0,'pm','l_0','мм');
printresult(fit2.x_0,fit2.d_x_0,'pm','x_0','мм');
printresult(fit2.y_0,fit2.d_y_0,'pm','y_0','мм');
printresult(fit2.gof.rmse,2,'s','RMSE','мВ');
%% ln(V/V_0) = f(x)
figure('Name','ln(V/V_0) = f(x)','Units','centimeters');
hold on;

X = [0 7];
plot(X,fit1.f(X),'Color','g');
plot(x,log(V./V(1)),'.','Color','b');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, мм','FontName','Liberation Serif');
ylabel('ln(\itV\rm/\itV\rm_0)','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% V(x)
figure('Name','V(x)','Units','centimeters');
hold on;

X = 0:0.1:7;
plot(X,fit2.f(X),'Color','g');
plot(x,V,'.','Color','b');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, мм','FontName','Liberation Serif');
ylabel('\itV\rm, мВ','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% Время жизни
img = imread('img\V(t).png');
t_per_px = 80./2314; % мкс/пикс
V_per_px = 600./1742; % мВ/пикс
fit3.range = [26 50]; % мкс

% Извлечение V(t) из изображения
t = t_per_px.*(1:size(img,2)).';
V = zeros(size(img,2),1);
for j = 1:size(img,2)
    V(j) = mean(find(img(:,j)));
end
V = V_per_px.*(max(V) - V);
t = t(~isnan(V));
V = V(~isnan(V));
% Аппроксимация ln(V) = f(t)
fit3.x = t(t > fit3.range(1) & t < fit3.range(2));
fit3.y = log(V(t > fit3.range(1) & t < fit3.range(2)));
fit3.f = fit(fit3.x,fit3.y,{'x','1'});
[fit3.a, fit3.b] = matsplit(coeffvalues(fit3.f));
[fit3.d_a, fit3.d_b] = matsplit(coeffvariances(fit3.f));
% τ, мкс
tau = -1./fit3.a;
d_tau = abs(fit3.d_a./fit3.a.^2);
%% Вывод результатов
printresult(fit3.a,fit3.d_a,'()','a','мкс^−1');
printresult(tau,d_tau,'()','τ','мкс');
%% V(t)
figure('Name','V(t)','Units','centimeters');
plot(t,V,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itt\rm, мкс','FontName','Liberation Serif');
ylabel('\itV\rm, мВ','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% ln(V) = f(t)
figure('Name','ln(V) = f(t)','Units','centimeters');
plot(fit3.x,fit3.y,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itt\rm, мкс','FontName','Liberation Serif');
ylabel('ln(\itV\rm), ln(мВ)','FontName','Liberation Serif');
xlim(fit3.range);
pbaspect([4 3 1]);
grid on;
figresize('width',17);