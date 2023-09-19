load('data_5_NI');
%% 1.1.
subpart = 1;
unpackstruct(part11(subpart));

Z_C.exp = U_C./(2.*sqrt(2).*I.*1e-3); % Ом
Z_C.theor = 1./(2.*pi.*f.*C.*1e-6); % Ом
tan_phi = Z_C.exp./R;
phi = atand(tan_phi); % °
if subpart == 1
    a = 1./(2.*pi.*f).*1e6; % Ом⋅мкФ
    % Аппроксимация Z_C(C^−1)
    fit1.f = fit(1./C,Z_C.exp,{'x'});
    fit1.a = coeffvalues(fit1.f);
    fit1.d_a = coeffvariances(fit1.f);
end
%% Вывод результатов
fprintf('\n\n1.1.%u.\n',subpart);
printresult(Z_C.exp,4,'s','Z_C (эксп.)','Ом');
printresult(Z_C.theor,4,'s','Z_C (теор.)','Ом');
printresult(tan_phi,4,'s','tan(φ)');
printresult(phi,4,'s','φ','°');
if subpart == 1
    printresult(a,4,'s','a (теор.)','Ом⋅мкФ');
    printresult(fit1.a,fit1.d_a,'pm','a (эксп.)','Ом⋅мкФ');
end
%% Z_C(C)
figure('Name','Z_C(C)','Units','centimeters');
hold on;

X = 0:0.002:1;
Y = a./X;
plot(X,Y,'Color','g');
plot(C,Z_C.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itC\rm, мкФ','FontName','Liberation Serif');
ylabel('\itZ_C\rm, Ом','FontName','Liberation Serif');
legend('Теоретическая кривая','Экспериментальные точки','Location','northeast');
ylim([0 3500]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% Z_C(C^−1)
figure('Name','Z_C(C^−1)','Units','centimeters');
hold on;

X = [0 100];
Y = a.*X;
plot(X,Y,'Color','g');
plot(1./C,Z_C.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(X,fit1.f(X),'Color','b');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itC\rm^{−1}, мкФ^{−1}','FontName','Liberation Serif');
ylabel('\itZ_C\rm, Ом','FontName','Liberation Serif');
legend('теоретическая прямая','экспериментальные точки','аппроксимирующая прямая','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 1.2.
subpart = [1,1];
unpackstruct(part12(subpart(1),subpart(2)));

switch subpart(1)
    case 1
        func_K = @(f) 1./sqrt(1 + (2.*pi.*f.*R.*C.*1e-6).^2);
        func_phi = @(f) -atand(2.*pi.*f.*R.*C.*1e-6);
    case 2
        func_K = @(f) 1./sqrt(1 + (2.*pi.*f.*R.*C.*1e-6).^-2);
        func_phi = @(f) atand(1./(2.*pi.*f.*R.*C.*1e-6));
end
% Теоретические зависимости
f.theor = logspace(2,4,500);
K.theor = func_K(f.theor);
L.theor = 20.*log10(K.theor);
phi.theor = func_phi(f.theor);
% Теоретические граничные значения
f_cutoff.theor = 1./(2.*pi.*R.*C.*1e-6);
K_cutoff.theor = func_K(f_cutoff.theor);
L_cutoff.theor = 20.*log10(K_cutoff.theor);
phi_cutoff.theor = func_phi(f_cutoff.theor);
% Экспериментальные граничные значения
points = K.exp > 0.1 & K.exp < 0.9;
f_cutoff.exp = interp1(K.exp(points),f.exp(points),1/sqrt(2));
K_cutoff.exp = interp1(f.exp,K.exp,f_cutoff.theor);
L_cutoff.exp = interp1(f.exp,L.exp,f_cutoff.theor);
phi_cutoff.exp = interp1(f.exp,phi.exp,f_cutoff.exp);
%% Вывод результатов
fprintf('\n\n1.2.%u.%u.\n',subpart);
fprintf('\nТеоретические значения\n');
printresult(f_cutoff.theor,4,'s','f_гр','Гц');
printresult(K_cutoff.theor,4,'s','K_гр');
printresult(L_cutoff.theor,4,'s','L_гр','дб');
printresult(phi_cutoff.theor,4,'s','φ_гр','°');
fprintf('\nЭкспериментальные значения\n');
printresult(f_cutoff.exp,4,'s','f_гр','Гц');
printresult(K_cutoff.exp,4,'s','K_гр');
printresult(L_cutoff.exp,4,'s','L_гр','дб');
printresult(phi_cutoff.exp,4,'s','φ_гр','°');
%% K(f)
figure('Name','K(f)','Units','centimeters');
hold on;

plot(f.theor,K.theor,'Color','g');
plot(f.exp,K.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'XScale','log','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(sprintf('\\itR\\rm\\bf = %g Ом, \\itC\\rm\\bf = %g мкФ',R,C),'FontName','Liberation Serif');
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itK','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','southwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% L(f)
figure('Name','L(f)','Units','centimeters');
hold on;

plot(f.theor,L.theor,'Color','g');
plot(f.exp,L.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'XScale','log','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(sprintf('\\itR\\rm\\bf = %g Ом, \\itC\\rm\\bf = %g мкФ',R,C),'FontName','Liberation Serif');
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itL\rm, дб','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','southwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% φ(f)
figure('Name','φ(f)','Units','centimeters');
hold on;

plot(f.theor,phi.theor,'Color','g');
plot(f.exp,phi.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'XScale','log','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(sprintf('\\itR\\rm\\bf = %g Ом, \\itC\\rm\\bf = %g мкФ',R,C),'FontName','Liberation Serif');
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itφ\rm, °','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','southwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2
subpart = 1;
unpackstruct(part2(subpart));

f_cutoff = (R + r)./(2.*pi.*L.*1e-3); % Гц
Z_L = U_L./(2.*sqrt(2).*I.*1e-3); % Ом
X_L.exp = sqrt(Z_L.^2 - r.^2); % Ом
X_L.theor = 2.*pi.*f.*L.*1e-3; % Ом
tan_phi = X_L.exp./(R + r);
phi = atand(tan_phi); % °
%% Вывод результатов
fprintf('\n\n2.%u.\n',subpart);
printresult(f_cutoff,4,'s','f_гр','Гц');
printresult(X_L.exp,4,'s','X_L (эксп.)','Ом');
printresult(X_L.theor,4,'s','X_L (теор.)','Ом');
printresult(tan_phi,4,'s','tan(φ)');
printresult(phi,4,'s','φ','°');
%% 3
unpackstruct(part3);

f_r = 1./(2.*pi.*sqrt(L.*C.*1e-9)); % Гц
Z_C = U_C./(2.*sqrt(2).*I.*1e-3); % Ом
Z_L = U_L./(2.*sqrt(2).*I.*1e-3); % Ом
X_L = sqrt(Z_L.^2 - r.^2); % Ом
tan_phi = (X_L - Z_C)./(R + r);
phi = atand(tan_phi); % °
%% Вывод результатов
fprintf('\n\n3.\n');
printresult(f_r,4,'s','f_р','Гц');
printresult(tan_phi,4,'s','tan(φ)');
printresult(phi,4,'s','φ','°');
%% U_R(f), U_C(f), U_L(f)
figure('Name','U_R(f), U_C(f), U_L(f)','Units','centimeters');
hold on;

X = logspace(log10(f(1)),log10(f(end)),500);
y = {U_R, U_C, U_L};
marker = {'o','s','^'};
color = {'b','r','g'};
label = {'U_R','U_C','U_L'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(X,pchip(f,y{j},X),'Color',color{j});
    subset(j) = plot(f,y{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',label{j});
end

set(gca,'XScale','log','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itU\rm, В','FontName','Liberation Serif');
legend(subset,'Location','northwest','FontName','Liberation Serif','FontAngle','italic');
pbaspect([4 3 1]);
grid on;
figresize('width',17);