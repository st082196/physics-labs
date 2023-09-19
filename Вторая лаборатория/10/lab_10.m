load('data_10');

% a_0, мкм
[a_0, d_a_0] = samplemean(sample_a_0,d_a_syst);
% x_0, мм
[x_0, d_x_0] = samplemean(sample_x_0,d_x_syst);
% Вычет нулевых отсчётов
a = abs(a - a_0);
x = abs(x - x_0);
sample_x_min = abs(sample_x_min - x_0);
% a′(a), мкм
a_img_a = abs(z_r_a - z_l_a);
% Аппроксимация a′(a)
fit1 = struct('f',cell(1,nnz(a<50)));
for i = 1:numel(fit1)
    fit1(i).f = fit(a(i:end),a_img_a(i:end),{'x','1'});
    [fit1(i).A, fit1(i).B] = matsplit(coeffvalues(fit1(i).f)); % ~, мкм
    [fit1(i).d_A, fit1(i).d_B] = matsplit(coeffvariances(fit1(i).f));
end
[~, i] = min([fit1.d_A]);
fit1 = fit1(i);
fit1.first_point = i;
% a_0 из аппроксимации, мкм
fit1.a_0 = -fit1.B./fit1.A;
fit1.d_a_0 = abs(fit1.a_0).*sqrt((fit1.d_A./fit1.A).^2 + (fit1.d_B./fit1.B).^2);
% a′(x), мкм
a_img_x = abs(z_r_x - z_l_x);
% Аппроксимация a′(x^−1)
fit2.x = x(1:fit2.num_points).^-1;
fit2.y = a_img_x(1:fit2.num_points);
fit2.f = fit(fit2.x,fit2.y,{'x'});
% x_min, мм
d_x_min_syst = sqrt(d_x_syst.^2 + d_x_0.^2);
[x_min, d_x_min] = samplemean(sample_x_min,d_x_min_syst);
% z_1, мкм
[z_1, d_z_1] = samplemean(sample_z_1,d_z_syst);
% z_2, мкм
[z_2, d_z_2] = samplemean(sample_z_2,d_z_syst);
% R_λ
R_wl = mean(wl)./diff(wl);
% R_x
R_x = abs(x_min./f_2.*(z_2-z_1)./diff(wl)).*1e4;
d_R_x = abs(R_x).*sqrt((d_x_min./x_min).^2 + (d_z_1./(z_2-z_1)).^2 + (d_z_2./(z_2-z_1)).^2);
%% Вывод результатов
printresult(a_0,d_a_0,'pm','a_0','мкм');
printresult(x_0,d_x_0,'pm','x_0','мм');
printresult(a,0,'d','a','мкм');
printresult(a_img_a,0,'d','a′','мкм');
printresult(fit1.A,fit1.d_A,'pm','A');
printresult(fit1.B,fit1.d_B,'pm','B','мкм');
printresult(fit1.a_0,fit1.d_a_0,'pm','a_0','мкм');
printresult(x,1,'d','x','мм');
printresult(a_img_x,0,'d','a′','мкм');
printresult(x_min,d_x_min,'pm','x','мм');
printresult(z_1,d_z_1,'pm','z_1','мкм');
printresult(z_2,d_z_2,'pm','z_2','мкм');
printresult(R_wl,0,'d','R_λ');
printresult(R_x,d_R_x,'pm','R_x');
%% a′(a)
figure('Name','a′(a)','Units','centimeters');
hold on;

X = [a(fit1.first_point) 300];
plot(X,fit1.f(X),'Color','g');
plot(a,a_img_a,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\ita\rm, мкм','FontName','Liberation Serif');
ylabel('\ita′\rm, мкм','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% a′(x)
figure('Name','a′(x)','Units','centimeters');
plot(x,a_img_x,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, мм','FontName','Liberation Serif');
ylabel('\ita′\rm, мкм','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('height',10);
%% a′(x^−1)
figure('Name','a′(x^−1)','Units','centimeters');
hold on;

X = [0 fit2.x(end) 6];
plot(X(1:2),[fit2.y(end) fit2.y(end)],'Color','g');
plot(X(2:3),fit2.f(X(2:3)),'Color','g');
plot(x.^-1,a_img_x,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm^{−1}, мм^{−1}','FontName','Liberation Serif');
ylabel('\ita′\rm, мкм','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);