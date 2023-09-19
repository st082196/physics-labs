load('data_128');

% R_x, Ом
R_x = mean(sample_R_x,2);
% d, мм
[d, d_d] = samplemean(sample_d,d_d_syst);
% Аппроксимация R_x(l)
fit_l_R_x = fit(l,R_x,{'x','1'});
[a, b] = matsplit(coeffvalues(fit_l_R_x)); % Ом/см, Ом
[d_a, d_b] = matsplit(coeffvariances(fit_l_R_x));
% ρ, Ом⋅м
rho = pi.*d.^2.*r./(4.*R_r).*a.*1e-4;
d_rho = abs(rho.*d_a./a);
%% Вывод результатов
printresult(d,d_d,'pm','d','мм');
printresult(a,d_a,'pm','a','Ом/см','auto');
printresult(b,d_b,'pm','b','Ом','auto');
printresult(rho,d_rho,'pm','ρ','Ом⋅м','auto');
%% R_x(l)
figure('Name','R_x(l)','Units','centimeters');
hold on;

X = [0 50];
plot(X,fit_l_R_x(X),'Color','g');
plot(l,R_x,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itl\rm, см','FontName','Liberation Serif');
ylabel('\itR_x\rm, Ом','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);