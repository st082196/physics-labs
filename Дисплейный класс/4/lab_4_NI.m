load('data_4_NI');
%% 1.
subpart = 1;
unpackstruct(part1(subpart));

I = U_pp/(2.*sqrt(2).*(R_1 + R_2)).*1e3; % мА
U_1 = U_pp.*R_1./(R_1 + R_2); % В
U_2 = U_pp.*R_2./(R_1 + R_2); % В
%% Вывод результатов
fprintf('\n\n1.%u.\n',subpart);
printresult(I,2,'d','I','мА');
printresult(U_1,3,'d','U_1','В');
printresult(U_2,3,'d','U_2','В');
%% 2.
subpart = 1;
unpackstruct(part2(subpart));

R_0 = R_1 + R_2.*R_3./(R_2 + R_3); % Ом
U_1 = U_pp.*R_1./R_0; % В
U_2_3 = U_pp.*(R_0 - R_1)./R_0; % В
I_1 = U_pp./(2.*sqrt(2).*R_0).*1e3; % мА
I_2 = U_pp.*(R_0 - R_1)./(2.*sqrt(2).*R_0.*R_2).*1e3; % мА
I_3 = U_pp.*(R_0 - R_1)./(2.*sqrt(2).*R_0.*R_3).*1e3; % мА
%% Вывод результатов
fprintf('\n\n2.%u.\n',subpart);
printresult(U_1,3,'d','U_1','В');
printresult(U_2_3,3,'d','U_2,3','В');
printresult(I_1,2,'d','I_1','мА');
printresult(I_2,2,'d','I_2','мА');
printresult(I_3,2,'d','I_3','мА');
%% 3.1.1.
unpackstruct(part311);

Z_C.exp = U_C./(2.*sqrt(2).*I_C.*1e-3); % Ом
Z_C.theor = 1./(2.*pi.*f.*C.*1e-6); % Ом
% Аппроксимация Z_C^−1(f)
fit1 = fit(f,1./Z_C.exp,{'x'});
a = coeffvalues(fit1);
d_a = coeffvariances(fit1);
% C (эксп), мкФ
C_exp = a./(2.*pi).*1e6;
d_C_exp = d_a./(2.*pi).*1e6;
%% Вывод результатов
fprintf('\n\n3.1.1.\n');
printresult(Z_C.exp,4,'s','Z_C, эксп','Ом');
printresult(Z_C.theor,4,'s','Z_C, теор','Ом');
printresult(1./Z_C.exp,4,'s','Z_C, эксп^−1','Ом^−1','auto');
printresult(a,d_a,'pm','a','Ом^−1⋅Гц^−1','auto');
printresult(C_exp,d_C_exp,'pm','C','мкФ');
%% Z_C(f)
figure('Name','Z_C(f)','Units','centimeters');
plot(f,Z_C.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itZ_C\rm, Ом','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% Z_C^−1(f)
figure('Name','Z_C^−1(f)','Units','centimeters');
hold on;

X = [0 5000];
plot(X,fit1(X),'Color','g');
plot(f,1./Z_C.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itZ_C\rm^{−1}, Ом^{−1}','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.1.2.
unpackstruct(part312);

Z_C.exp = U_C./(2.*sqrt(2).*I_C.*1e-3); % Ом
Z_C.theor = 1./(2.*pi.*f.*C.*1e-6); % Ом
%% Вывод результатов
fprintf('\n\n3.1.2.\n');
for j = 1:numel(source)
    fprintf('\n%s\n',source(j));
    printresult(Z_C.exp(:,j),4,'s','Z_C, эксп','Ом');
    printresult(Z_C.theor(:,j),4,'s','Z_C, теор','Ом');
end
%% Z_C(С)
j = 1;

figure('Name','Z_C(С)','Units','centimeters');
plot(C,Z_C.exp(:,j),'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
title(source(j));
xlabel('\itC\rm, мкФ','FontName','Liberation Serif');
ylabel('\itZ_C\rm, Ом','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.2.
unpackstruct(part32);

Z_C_1 = 1./(2.*pi.*f.*C_1); % Ом
Z_C_2 = 1./(2.*pi.*f.*C_2); % Ом
I = U_pp./(2.*sqrt(2).*(Z_C_1 + Z_C_2)).*1e3; % мА
U_1 = U_pp.*Z_C_1./(Z_C_1 + Z_C_2); % В
U_2 = U_pp.*Z_C_2./(Z_C_1 + Z_C_2); % В
%% Вывод результатов
fprintf('\n\n3.2.\n');
printresult(I,2,'d','I','мА');
printresult(U_1,3,'d','U_1','В');
printresult(U_2,3,'d','U_2','В');
%% I(f)
figure('Name','I(f)','Units','centimeters');
hold on;

X = [0 5000];
Y = pi.*X.*U_pp./(sqrt(2).*(1./C_1 + 1./C_2)).*1e3;
plot(X,Y,'Color','g');
plot(f,I,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itI\rm, мА','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 4.
subpart = 1;
unpackstruct(part4(subpart));

Y_C_1 = 2.*pi.*f.*C_1; % Ом
Y_C_2 = 2.*pi.*f.*C_2; % Ом
I = U_pp.*(Y_C_1 + Y_C_2)./(2.*sqrt(2)).*1e3; % мА
I_1 = U_pp.*Y_C_1./(2.*sqrt(2)).*1e3; % мА
I_2 = U_pp.*Y_C_2./(2.*sqrt(2)).*1e3; % мА
%% Вывод результатов
fprintf('\n\n4.%u.\n',subpart);
printresult(I,2,'d','I','мА');
printresult(I_1,2,'d','I_1','мА');
printresult(I_2,2,'d','I_2','мА');
%% 5.
subpart = 1;
unpackstruct(part5(subpart));


L = D.^2.*N.^2./(45.*D + 100.*l).*1e-3; % мГн
I = U_R./R; % А
Z = U_in./(2.*sqrt(2).*I); % Ом
Z_L = U_L./I; % Ом
% Аппроксимация Z^2(f^2)
[x, y, fit2] = deal(cell(size(R)));
[a, b, d_a, d_b] = deal(zeros(size(fit2)));
for j = 1:numel(R)
    x{j} = f(1:num_fit_points(j),j).^2;
    y{j} = Z(1:num_fit_points(j),j).^2;
    fit2{j} = fit(x{j},y{j},{'x','1'});
    [a(j), b(j)] = matsplit(coeffvalues(fit2{j}));
    [d_a(j), d_b(j)] = matsplit(coeffvariances(fit2{j}));
end
% L (эксп), мГн
L_exp = sqrt(a)./(2.*pi).*1e3;
d_L_exp = d_a./(4.*pi.*sqrt(a)).*1e3;
% r, Ом
r = sqrt(b) - R;
d_r = d_b./(2.*sqrt(b));
U_X_L = I.*2.*pi.*f.*L_exp.*1e-3; % В
%% Вывод результатов
fprintf('\n\n5.%u.\n',subpart);
printresult(L,1,'d','L','мГн');
for j = 1:numel(R)
    fprintf('\nR = %u Ом\n',R(j));
    printresult(Z(:,j),4,'s','Z','Ом');
    printresult(Z_L(:,j),4,'s','Z_L','Ом');
    printresult(U_X_L(:,j),3,'d','U_X_L','В');
    printresult(a(j),d_a(j),'pm','a','Ом^2⋅Гц^−2','auto');
    printresult(b(j),d_b(j),'pm','b','Ом^2');
    printresult(L_exp(j),d_L_exp(j),'pm','L','мГн');
    printresult(r(j),d_r(j),'pm','r','Ом');
end
%% Z^2(f^2)
j = 1;

figure('Name','Z^2(f^2)','Units','centimeters');
hold on;

X = [0 2.5e6];
plot(X,fit2{j}(X),'Color','g');
plot(x{j},y{j},'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(sprintf('\\itN\\rm\\bf = %u, \\itR\\rm\\bf = %g Ом',N,R(j)),'FontName','Liberation Serif');
xlabel('\itf\rm^2, Гц^2','FontName','Liberation Serif');
ylabel('\itZ\rm^2, Ом^2','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);