load('data_14');

% τ, с
[tau, d_tau] = samplemean(sample_tau);
% K, мПа⋅см^3⋅г^−1
K = eta_ref(1)./((rho_b-rho_l(1)).*tau(1));
d_K = abs(K.*d_tau(1)./tau(1));
% η, мПа⋅с
eta = [eta_ref(1); K.*(rho_b-rho_l(2:end)).*tau(2:end)];
d_eta = [0.0001; abs(eta(2:end)).*sqrt((d_K./K).^2 + (d_tau(2:end)./tau(2:end)).^2)];
% Аппроксимация y(x)
x = 1./R.*(1./(T+273.15) - 1./(T(1)+273.15)); % моль/Дж
y = log(eta./eta(1));
fit_x_y = fit(x,y,{'x'});
% W, Дж/моль
W = coeffvalues(fit_x_y);
d_W = coeffvariances(fit_x_y);
%% Вывод результатов
printresult(tau,d_tau,'pm','τ','с');
printresult(K,d_K,'pm','K','мПа⋅см^3⋅г^−1','auto');
printresult(eta,d_eta,'pm','η','мПа⋅с');
printresult(W,d_W,'pm','W','Дж/моль');
%% η(T)
figure('Name','η(T)','Units','centimeters');
hold on;

errorbar(T,eta,d_eta,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'DisplayName','экспериментальные значения');
plot(T,eta_ref,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'DisplayName','табличные значения');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itT\rm, °C','FontName','Liberation Serif');
ylabel('\itη\rm, мПа \cdot с','FontName','Liberation Serif');
legend('Location','southwest');
pbaspect([4 3 1]);
grid on;
figresize('width',10.5);
%% y(x)
figure('Name','y(x)','Units','centimeters');
hold on;

X = [-3.5e-5 0];
plot(X,fit_x_y(X),'Color','g');
plot(x,y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, моль/Дж','FontName','Liberation Serif');
ylabel('\ity','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);