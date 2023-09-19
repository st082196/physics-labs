load('data_21');
%% 1.
unpackstruct(part1);

% Аппроксимация l(n)
fit_n_l = fit(n,l,{'x','1'});
% λ, мкм
wl = matsplit(coeffvalues(fit_n_l)).*1e4;
d_wl = matsplit(coeffvariances(fit_n_l)).*1e4;
% c, м/с
c = wl.*f.*1e-3;
d_c = sqrt((f.*d_wl).^2 + (wl.*d_f).^2).*1e-3;
%% Вывод результатов
printresult(wl,d_wl,'pm','λ','мкм');
printresult(c,d_c,'pm','c','м/с');
%% l(n)
figure('Name','l(n)','Units','centimeters');
hold on;

x_fit = [0 200];
plot(x_fit,fit_n_l(x_fit),'Color','g');
plot(n,l,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itn','FontName','Liberation Serif');
ylabel('\itl\rm, см','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);
i = [1:numel(l)-1; 2:numel(l)].';
% i = [1, 4];

c = abs((l(i(:,2)) - l(i(:,1)))./(m(i(:,2))./(f_m(i(:,2)) - f(i(:,2))) - m(i(:,1))./(f_m(i(:,1)) - f(i(:,1))))).*10;
d_c = abs(c).*sqrt( ...
    2.*(d_l./(l(i(:,2)) - l(i(:,1)))).^2 + ...
    2.*( ...
    (m(i(:,1)).*(f_m(i(:,2)) - f(i(:,2)))./(f_m(i(:,1)) - f(i(:,1)))).^2 + ...
    (m(i(:,2)).*(f_m(i(:,1)) - f(i(:,1)))./(f_m(i(:,2)) - f(i(:,2)))).^2 ...
    ).*(d_f./(m(i(:,2)).*(f_m(i(:,1)) - f(i(:,1))) - m(i(:,1)).*(f_m(i(:,2)) - f(i(:,2))))).^2 ...
    );
mean_c = mean(c);
d_mean_c = sqrt(sum(d_c.^2))./numel(c);
%% Вывод результатов
printresult(c,d_c,'pm','c','м/с');
printresult(mean_c,d_mean_c,'pm','mean_c','м/с');
%% 3.
unpackstruct(part3);

% Аппроксимация l(t)
fit_t_l = fit(t,l,{'x','1'});
% u, м/с
u = abs(matsplit(coeffvalues(fit_t_l))).*1e4;
d_u = abs(matsplit(coeffvariances(fit_t_l))).*1e4;
%% Вывод результатов
printresult(u,d_u,'pm','u','м/с');
%% l(t)
figure('Name','l(t)','Units','centimeters');
hold on;

x_fit = [10 80];
plot(x_fit,fit_t_l(x_fit),'Color','g');
plot(t,l,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itt\rm, мкс','FontName','Liberation Serif');
ylabel('\itl\rm, см','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 4.
unpackstruct(part4);

% λ, мкм
wl = m.*WL.*L./d.*1e-3;
d_wl = abs(wl).*sqrt((d_L./L).^2 + (d_d./d).^2);
%% Вывод результатов
printresult(wl,d_wl,'pm','λ','мкм');