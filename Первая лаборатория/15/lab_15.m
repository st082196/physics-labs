load('data_15');
%% 1. & 2.
unpackstruct(part1);

% F, мН
[F, d_F] = samplemean(sample_F,d_F_syst);
% σ, мН/м
sigma = F./(pi.*(D_1 + D_2)).*1e3;
d_sigma = abs(sigma).*sqrt(2.*(d_D./(D_1+D_2)).^2 + (d_F./F).^2);
%% Вывод результатов
printresult(F,d_F,'pm','F','мН');
printresult(sigma,d_sigma,'pm','σ','мН/м');
%% σ(C)
figure('Name','σ(C)','Units','centimeters');
hold on;

X = 0:100;
plot(X,spline(C_ref,sigma_ref,X),'Color','g');
plot(C_ref,sigma_ref,'s','Color','g','MarkerFaceColor','g','MarkerSize',3);
plot(C,sigma,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('\itT\rm\bf = 20°C','FontName','Liberation Serif');
xlabel('\itC\rm, %','FontName','Liberation Serif');
ylabel('\itσ\rm, мН/м','FontName','Liberation Serif');
legend('','табличные значения','экспериментальные значения','Location','northeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.
unpackstruct(part3);

% σ, мН/м
sigma = F./(pi.*(D_1 + D_2)).*1e3;
d_sigma = abs(sigma).*sqrt(2.*(d_D./(D_1+D_2)).^2 + (d_F./F).^2);
% Аппроксимация σ(T)
fit1 = fit(T_ref,sigma_ref,{'x','1'});
fit2 = fit(T,sigma,{'x','1'});
[a, b] = matsplit(coeffvalues(fit2)); % мН/(м⋅°C), мН/м
[d_a, d_b] = matsplit(coeffvariances(fit2));
% K_0, Дж/(моль^(2/3)⋅К)
K_0 = -V_0.^(2./3).*a.*1e-5;
d_K_0 = abs(K_0.*d_a./a);
%% Вывод результатов
printresult(sigma,d_sigma,'pm','σ','мН/м');
printresult(a,d_a,'pm','a','мН/(м⋅°C)');
printresult(b,d_b,'pm','b','мН/м');
printresult(K_0,d_K_0,'pm','K_0','Дж/(моль^(2/3)⋅К)','auto');
%% σ(T)
figure('Name','σ(T)','Units','centimeters');
hold on;

X = [0 100];
plot(X,fit1(X),'Color','g');
plot(T_ref,sigma_ref,'s','Color','g','MarkerFaceColor','g','MarkerSize',3');
plot(X,fit2(X),'Color','b');
errorbar(T,sigma,d_sigma,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itT\rm, °C','FontName','Liberation Serif');
ylabel('\itσ\rm, мН/м','FontName','Liberation Serif');
legend('','табличные значения','','экспериментальные значения','Location','northeast');
pbaspect([4 3 1]);
grid on;
figresize('width',14.5);