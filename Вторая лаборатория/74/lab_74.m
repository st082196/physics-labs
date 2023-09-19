load('data_74');
part = 1; % Часть работы (1 или 2)

switch part
    case 1
        nu_12 = nu_12_V;
        T = T(1);
    case 2
        nu_12 = nu_12_T;
        V_2 = V_2(end);
end

% P_1, см^3/моль
P_1 = (eps_1 - 1)./(eps_1 + 2).*M_1./rho_1;
% ε_12
eps_12 = (2.*(nu_0-nu_12)./nu_12+((nu_0-nu_12)./nu_12).^2) .* (eps_1-1) .* (2.*(nu_0-nu_1)./nu_1+((nu_0-nu_1)./nu_1).^2)^-1 + 1;
% M_12, г/моль
M_12 = M_1.*M_2.*(rho_1.*V_1 + rho_2.*V_2)./(M_1.*rho_2.*V_2 + M_2.*rho_1.*V_1);
% ρ_12, г/см^3
rho_12 = (rho_1.*V_1 + rho_2.*V_2)./(V_1 + V_2);
% P_12, см^3/моль
P_12 = (eps_12 - 1)./(eps_12 + 2).*M_12./rho_12;
% P_2, см^3/моль
P_2 = (P_12.*(M_1.*rho_2.*V_2 + M_2.*rho_1.*V_1) - P_1.*M_2.*rho_1.*V_1)./(M_1.*rho_2.*V_2);
% R, см^3/моль
R = (n_2.^2 - 1)./(n_2.^2 + 2).*M_2./rho_2;
switch part
    case 1
        % μ, Д
        mu = 3./2.*sqrt((P_2-R).*k.*(T+273.15)./(pi.*N_A)).*1e18;
        [mean_mu, d_mu] = samplemean(mu);
    case 2
        % Аппроксимация P_2 − R = f(T^−1)
        x = (T+273.15).^-1;
        y = P_2 - R;
        fit_x_y = fit(x,y,{'x','1'});
        [a, b] = matsplit(coeffvalues(fit_x_y)); % см^3⋅К⋅моль^−1, см^3/моль
        [d_a, d_b] = matsplit(coeffvariances(fit_x_y));
        % μ, Д
        mu = 3./2.*sqrt(a.*k./(pi.*N_A)).*1e18;
        d_mu = abs(mu./2.*d_a./a);
end
%% Вывод результатов
printresult(P_1,3,'s','P_1','см^3/моль');
printresult(eps_12,3,'s','ε_12');
printresult(M_12,4,'s','M_12','г/моль');
printresult(rho_12,3,'s','ρ_12','г/см^3');
printresult(P_12,3,'s','P_12','см^3/моль');
printresult(P_2,3,'s','P_2','см^3/моль');
printresult(R,3,'s','R','см^3/моль');
switch part
    case 1
        printresult(mu,3,'s','μ','Д');
        printresult(mean_mu,d_mu,'pm','mean_μ','Д');
    case 2
        printresult(y,3,'s','P_2 − R','см^3/моль');
        printresult(a,d_a,'pm','a','см^3⋅К⋅моль^−1');
        printresult(mu,d_mu,'pm','mean_μ','Д');
end
%% P_2 − R = f(T^−1)
figure('Name','P_2 − R = f(T^−1)','Units','centimeters');
hold on;

X = [3 3.5].*1e-3;
plot(X,fit_x_y(X),'Color','g');
plot(x,y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itT\rm^{−1}, K^{−1}','FontName','Liberation Serif');
ylabel('\itP\rm_2 − \itR\rm, см^3/моль','FontName','Liberation Serif');
xlim(X);
pbaspect([4 3 1]);
grid on;
figresize('width',17);