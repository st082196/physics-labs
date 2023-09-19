load('data_4');

[theta, Q, x, y, fit_x_y, c, d_c, V_inv_R, V_inv_rho] = deal(cell(size(doubletheta)));
for i = 1:numel(doubletheta)
    theta{i} = doubletheta{i}./2;
    Q{i} = sind(theta{i}).^2./sind(theta{i}(1)).^2;
    % Аппроксимация y(x)
    x{i} = sqrt(x_sqr{i});
    y{i} = sind(theta{i});
    fit_x_y{i} = fit(x{i},y{i},{'x'});
    c{i} = coeffvalues(fit_x_y{i});
    d_c{i} = coeffvariances(fit_x_y{i});
    % a, Å
    a.(name{i}) = wl./(2.*c{i});
    d_a.(name{i}) = abs(a.(name{i}).*d_c{i}./c{i});
end
% R, Å
R.Br = a.KBr./2 - R.K;
R.I = sqrt(3)./2.*a.KI - R.K;
R.Cs(1) = sqrt(3)./2.*a.CsCl - R.Cl;
R.Cs(2) = sqrt(3)./2.*a.CsBr - R.Br;
R.Cs(3) = sqrt(3)./2.*a.CsI - R.I;
% V^−1(R_+,R_−)
V_inv_R{1} = 3.*sqrt(3)./(8.*(R.Cs(1) + R.Cl).^3);
V_inv_R{2} = 3.*sqrt(3)./(8.*(R.Cs(2) + R.Br).^3);
V_inv_R{3} = 3.*sqrt(3)./(8.*(R.Cs(3) + R.I).^3);
V_inv_R{4} = 1./(2.*(R.K + R.Cl).^3);
V_inv_R{5} = 1./(2.*(R.K + R.Br).^3);
V_inv_R{6} = 3.*sqrt(3)./(8.*(R.K + R.I).^3);
% V^−1(ρ), Å^−3
V_inv_rho{1} = rho{1}./((m.Cs + m.Cl).*Da).*1e-27;
V_inv_rho{2} = rho{2}./((m.Cs + m.Br).*Da).*1e-27;
V_inv_rho{3} = rho{3}./((m.Cs + m.I).*Da).*1e-27;
V_inv_rho{4} = rho{4}./((m.K + m.Cl).*Da).*1e-27;
V_inv_rho{5} = rho{5}./((m.K + m.Br).*Da).*1e-27;
V_inv_rho{6} = rho{6}./((m.K + m.I).*Da).*1e-27;
%% Вывод результатов
for i = 1:numel(doubletheta)
    fprintf('\n%s\n',name{i});
    printresult(doubletheta{i},2,'d','2θ','°');
    printresult(Q{i},2,'d','Q');
    printresult(c{i},d_c{i},'()','c');
    printresult(a.(name{i}),d_a.(name{i}),'pm','a','Å');
end
printresult(cell2mat(V_inv_R),3,'s','V^−1(R_+,R_−)','Å^−3');
printresult(cell2mat(V_inv_rho),3,'s','V^−1(ρ)','Å^−3');
printresult(R.Br,3,'s','R_Br','Å');
printresult(R.I,3,'s','R_I','Å');
printresult(R.Cs,3,'s','R_Cs','Å');
%% y(x)
figure('Name','y(x)','Units','centimeters');
hold on;

marker = {'o','s','^','v','p','x'};
color = {'b','r','g','m','k','#0072BD'};
subset = zeros(size(color));
for i = 1:numel(subset)
    plot(x{i},fit_x_y{i}(x{i}),'Color',color{i});
    subset(i) = plot(x{i},y{i},marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('x','FontName','Liberation Serif','FontAngle','italic');
ylabel('y','FontName','Liberation Serif','FontAngle','italic');
legend(subset,name,'Location','southeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);