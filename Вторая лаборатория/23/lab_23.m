load('data_23');

% ν, см^−1
nu = 1./(n.*wl.*1e-7);
% Аппроксимация ν(m) по одиночным линиям
fit_m_nu_half = fit(m(1:num_single_lines),nu(i_m(1:num_single_lines)),'poly2');
p_half = coeffvalues(fit_m_nu_half); % см^−1
d_p_half = coeffvariances(fit_m_nu_half);
% Аппроксимация ν(m) по всем линиям
fit_m_nu = fit(m,nu(i_m),'poly2');
p = coeffvalues(fit_m_nu); % см^−1
d_p = coeffvariances(fit_m_nu);
% B, см^−1
B(1) = (p(1) + p(2))./2;
d_B(1) = sqrt(d_p(1).^2 + d_p(2).^2)./2;
B(2) = (p(2) - p(1))./2;
d_B(2) = sqrt(d_p(1).^2 + d_p(2).^2)./2;
% I, кг⋅м^2
I = h./(8.*pi.^2.*B.*1e2.*c);
d_I = abs(I.*d_B./B);
% μ, а. е. м.
mu = m_C.*m_N./(m_C + m_N);
% r, Å
r = sqrt(I./(mu.*Da)).*1e10;
d_r = abs(d_I./(2.*sqrt(mu.*Da.*I))).*1e10;
%% Вывод результатов
printresult(wl(i_m),6,'s','λ','нм');
printresult(nu(i_m),6,'s','ν','см^−1');
printresult(p_half,d_p_half,'()','p, аппроксимация по одиночным линиям','см^−1');
printresult(p,d_p,'()','p, аппроксимация по всем линиям','см^−1');
printresult(B,d_B,'()','B','см^−1');
printresult(I,d_I,'()','I','кг⋅м^2','auto');
printresult(r,d_r,'()','r','Å');
%% I(λ)
figure('Name','I(λ)','Units','centimeters');
plot(wl,intensity,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
xlim([387 388.4]);
ylim([0 1200]);
pbaspect([2 1 1]);
figresize('width',16);
%% I(ν)
figure('Name','I(ν)','Units','centimeters');
plotobject = plot(nu,intensity,'Color','b');
plotobject.DataTipTemplate.DataTipRows(3).Label = 'Index';
plotobject.DataTipTemplate.DataTipRows(3).Value = 1:numel(plotobject.XData);
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itν\rm, см^{−1}','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
xlim([25740 25830]);
ylim([0 1200]);
pbaspect([2 1 1]);
figresize('width',16);
%% ν(m), аппроксимация по одиночным линиям
figure('Name','ν(m), аппроксимация по одиночным линиям','Units','centimeters');
hold on;

X = -62:5;
plot(X,fit_m_nu_half(X),'Color','g');
plot(m(1:num_single_lines),nu(i_m(1:num_single_lines)),'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itm','FontName','Liberation Serif');
ylabel('\itν\rm, см^{−1}','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% ν(m), аппроксимация по всем линиям
figure('Name','ν(m), аппроксимация по всем линиям','Units','centimeters');
hold on;

plot(-62:5,fit_m_nu(-62:5),'Color','g');
plot(m,nu(i_m),'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itm','FontName','Liberation Serif');
ylabel('\itν\rm, см^{−1}','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);