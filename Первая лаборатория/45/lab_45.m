load('data_45');
%% 1.
unpackstruct(part1);

% E_s, клк
[E_s, d_E_s] = samplemean(sample_E_s,d_E_s_syst);
%% Вывод результатов
printresult(E_s,d_E_s,'pm','E_s','клк');
%% 2.
unpackstruct(part2(1));

% Аппроксимация E_s(x)
fit_x_E_s = fit(x,E_s,'I_s/(x-x_0)^2','StartPoint',[1000,100],'Robust','Bisquare');
[I_s, x_0] = matsplit(coeffvalues(fit_x_E_s)); % клк⋅см^2, см
[d_I_s, d_x_0] = matsplit(coeffvariances(fit_x_E_s));
r = abs(x - x_0); % см
%% Вывод результатов
fprintf('\n%s\n',label);
printresult(I_s,d_I_s,'pm','I_s','клк⋅см^2');
printresult(x_0,d_x_0,'pm','x_0','см');
%% E_s(x)
figure('Name','E_s(x)','Units','centimeters');
hold on;

x_fit = 15:0.1:45;
plot(x_fit,fit_x_E_s(x_fit),'Color','g');
plot(x,E_s,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(label);
xlabel('\itx\rm, см','FontName','Liberation Serif');
ylabel('\itE\rm_s, клк','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% E_s(r^−2)
figure('Name','E_s(r^−2)','Units','centimeters');
hold on;

x_fit = [0 3].*1e-3;
plot(x_fit,I_s.*x_fit,'Color','g');
plot(r.^-2,E_s,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(label);
xlabel('\itr\rm^{−2}, см^{−2}','FontName','Liberation Serif');
ylabel('\itE\rm_s, клк','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.
unpackstruct(part3)

% E_s, клк
[E_s, d_E_s] = samplemean(sample_E_s,d_E_s_syst);
% Φ_s, лм
PHI_s = E_s.*pi.*d.^2./4.*1e-1;
d_PHI_s = abs(pi.*d./2).*sqrt((d./2.*d_E_s).^2 + (E_s.*d_d).^2).*1e-1;
% Φ_e, мВт
PHI_e = a.*(e - e_0);
d_PHI_e = abs(a).*sqrt(2.*d_e.^2);
% K
K = PHI_s./(V.*PHI_e).*1e3;
d_K = abs(K).*sqrt((d_PHI_s./PHI_s).^2 + (d_PHI_e./PHI_e).^2);
%% Вывод результатов
printresult(E_s,d_E_s,'pm','E_s','клк');
printresult(PHI_s,d_PHI_s,'pm','Φ_s','лм');
printresult(PHI_e,d_PHI_e,'pm','Φ_e','мВт');
printresult(K,d_K,'pm','K');