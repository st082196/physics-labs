load('data_33');

% α, см
alpha_I = (alpha_1_I + alpha_2_I)./2;
alpha_T = (alpha_1_T + alpha_2_T)./2;
alpha_0 = min(alpha_T);
% H, А/м
H = I_working.*N./l;
d_H = abs(N./l.*d_I);
% m
m = (1./eps.^2 - 1).*(1./(2.*eps).*log((1+eps)./(1-eps)) - 1);
% C′_b, Вб/м
C_b = 2.*mu_0.*H.*S_0.*N_0./alpha_0.*1e-2;
d_C_b = abs(C_b).*sqrt((d_H./H).^2 + (d_alpha./alpha_0).^2);
% B, Тл
B = C_b.*alpha_T./(2.*N_0.*S_Ni).*1e2 - mu_0.*H.*(S_0-S_Ni)./S_Ni;
d_B = sqrt( ...
    (alpha_T./(2.*N_0.*S_Ni).*d_C_b.*1e2).^2 + ...
    (C_b./(2.*N_0.*S_Ni).*d_alpha.*1e2).^2 + ...
    (mu_0.*(S_0-S_Ni)./S_Ni.*d_H).^2 ...
    );
% H′, А/м
H_corrected = (H - m.*B./mu_0)./(1 - m);
d_H_corrected = abs(1./(1 - m)).*sqrt(d_H.^2 + (m./mu_0.*d_B).^2);
% M, А/м
M = B./mu_0 - H_corrected;
d_M = sqrt((d_B./mu_0).^2 + d_H_corrected.^2);
% n, м^−3
n = N_A.*rho_Ni./m_Ni.*1e6;
% μ, А⋅м^2
mu = M(1)./n;
d_mu = abs(d_M(1)./n);
% γ, Дж⋅А^−2⋅м^−1
gamma = 3.*k.*(T_c+273.15)./(n.*mu.^2);
d_gamma = abs(gamma).*sqrt((d_T_c./(T_c+273.15)).^2 + (2.*d_mu./mu).^2);
%% Вывод результатов
printresult(alpha_I,2,'d','α(I)','см');
printresult(H,d_H,'pm','H','А/м');
printresult(alpha_T,2,'d','α(T)','см');
printresult(m,4,'s','m');
printresult(C_b,d_C_b,'pm','C′_b','Вб/м');
printresult(B,d_B,'pm','B','Тл');
printresult(H_corrected,d_H_corrected,'pm','H′','А/м',3);
printresult(M,d_M,'pm','M','А/м',3);
printresult(n,3,'s','n','м^−3','auto');
printresult(mu,d_mu,'pm','μ','А⋅м^2','auto');
printresult(gamma,d_gamma,'pm','γ','Дж⋅А^−2⋅м^−1','auto');
%% α(I)
figure('Name','α(I)','Units','centimeters');
plot(I,alpha_I,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itI\rm, А','FontName','Liberation Serif');
ylabel('\itα\rm, см','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% α(T)
figure('Name','α(T)','Units','centimeters');
plot(T,alpha_T,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itT\rm, °C','FontName','Liberation Serif');
ylabel('\itα\rm, см','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% M(T)
figure('Name','M(T)','Units','centimeters');
plot(T,M,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itT\rm, °C','FontName','Liberation Serif');
ylabel('\itM\rm, А/м','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);