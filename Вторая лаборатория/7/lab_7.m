load('data_7');

% λ, Å
wl = spline(phi_cal,wl_cal,phi);
% Φ, мкВт
PHI = U./eps;
d_PHI = abs(d_U./eps);
% K, А/Вт
K = I./PHI;
d_K = abs(K).*sqrt((d_I./I).^2 + (d_PHI./PHI).^2);
%% Вывод результатов
printresult(wl,4,'s','λ','Å');
printresult(U,d_U,'pm','U','мкВ');
printresult(PHI,d_PHI,'pm','Φ','мкВт');
printresult(I,d_I,'pm','I','мкА');
printresult(K,d_K,'pm','K','А/Вт');
%% λ(φ)
figure('Name','λ(φ)','Units','centimeters');
hold on;

X = min(phi_cal):10:max(phi_cal);
plot(X,spline(phi_cal,wl_cal,X),'Color','g');
plot(phi_cal,wl_cal,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itφ\rm, °','FontName','Liberation Serif');
ylabel('\itλ\rm, Å','FontName','Liberation Serif');
pbaspect([1 1 1]);
grid on;
figresize('width',10.5);
%% Φ(λ), I(λ), K(λ)
figure('Name','Φ(λ), I(λ), K(λ)','Units','centimeters');
hold on;

plot(wl,PHI./max(PHI),'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'DisplayName','Φ');
plot(wl,I./max(I),'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'DisplayName','I');
plot(wl,K./max(K),'^','Color','g','MarkerFaceColor','g','MarkerSize',3,'DisplayName','K');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, Å','FontName','Liberation Serif');
legend('Location','northwest','FontName','Liberation Serif','FontAngle','italic');
pbaspect([4 3 1]);
grid on;
figresize('width',17);