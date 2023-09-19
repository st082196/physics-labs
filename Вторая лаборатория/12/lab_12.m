load('data_12');

% x, мм
x = x_0.*N_spectrum.*1e-3;
% ∆λ, нм
delta_wl = diff(wl_spectrum);
% ∆x, мм
delta_x = diff(x);
% ⟨λ⟩, нм
wl_mean = (wl_spectrum(1:end-1) + wl_spectrum(2:end))./2;
% Λ, нм/мм
LAMBDA = delta_wl./delta_x;
% Аппроксимация λ(N)
fit_N_wl = polyfit(N_spectrum,wl_spectrum,3);
% Аппроксимация Λ(⟨λ⟩)
fit_wl_LAMBDA = polyfit(wl_mean,LAMBDA,2);
% λ, нм
wl = polyval(fit_N_wl,N);
% ν, см^−1
nu = 1e7./(n.*wl);
%% Вывод результатов
printresult(x,3,'d','x','мм');
printresult(delta_wl,3,'d','∆λ','нм');
printresult(delta_x,3,'d','∆x','мм');
printresult(wl_mean,3,'d','⟨λ⟩','нм');
printresult(LAMBDA,3,'d','Λ','нм/мм');
printresult(fit_N_wl,4,'s','Коэффициенты λ(N)',[],'auto-unique');
printresult(fit_wl_LAMBDA,4,'s','Коэффициенты Λ(⟨λ⟩)',[],'auto-unique');
%% I(N)
i = 3;

figure('Name',"I(N), d = "+d(i)+" мм",'Units','centimeters');
plot(N,I(:,i),'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title("\itd\rm\bf = "+d(i)+" мм",'FontName','Liberation Serif');
xlabel('\itN','FontName','Liberation Serif');
xlim([N(1) N(end)]);
ylim([0 4000]);
pbaspect([2 1 1]);
figresize('width',17);
%% λ(N)
figure('Name','λ(N)','Units','centimeters');
hold on;

X = N(1):10:N(end);
plot(X,polyval(fit_N_wl,X),'Color','g');
plot(N_spectrum,wl_spectrum,'.','Color','b');

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itN','FontName','Liberation Serif');
ylabel('\itλ\rm, нм','FontName','Liberation Serif');
xlim([N(1) N(end)]);
ylim([550 max(wl_spectrum)]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% Λ(⟨λ⟩)
figure('Name','Λ(⟨λ⟩)','Units','centimeters');
hold on;

X = wl(1):1:wl(end);
plot(X,polyval(fit_wl_LAMBDA,X),'Color','g');
plot(wl_mean,LAMBDA,'.','Color','b');

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\langle\itλ\rm\rangle, нм','FontName','Liberation Serif');
ylabel('\itΛ\rm, нм/мм','FontName','Liberation Serif');
xlim([wl(1) wl(end)]);
ylim([5 max(LAMBDA)]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% I(λ)
i = 1;

figure('Name',"I(λ), d = "+d(i)+" мм",'Units','centimeters');
plot(wl,I(:,i),'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title("\itd\rm\bf = "+d(i)+" мм",'FontName','Liberation Serif');
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
xlim([wl(1) wl(end)]);
ylim([0 4000]);
pbaspect([2 1 1]);
figresize('width',14);
%% I(ν)
i = 1;

figure('Name',"I(ν), d = "+d(i)+" мм",'Units','centimeters');
plot(nu,I(:,i),'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title("\itd\rm\bf = "+d(i)+" мм",'FontName','Liberation Serif');
xlabel('\itν\rm, см^{−1}','FontName','Liberation Serif');
xlim([nu(end) nu(1)]);
ylim([0 4000]);
pbaspect([2 1 1]);
figresize('width',14);