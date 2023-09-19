load('data_2');

% λ_0, пм
wl_0 = 2.*d.*sind(theta_0);
d_wl_0 = abs(2.*d.*cosd(theta_0).*deg2rad(d_theta_0));
% Аппроксимация λ_0(U^−1)
U_inv = U.^-1;
fit1 = fit(U_inv,wl_0,{'x'});
a = coeffvalues(fit1); % пм⋅кВ
d_a = coeffvariances(fit1);
% a (теор.), пм⋅кВ
a_theor = h.*c./e.*1e9;
% λ (L-серия), пм
LiF.wl_L = 2.*LiF.d.*sind(LiF.theta_L)./LiF.n_L;
LiF.d_wl_L = abs(2.*LiF.d.*cosd(LiF.theta_L)./LiF.n_L.*deg2rad(LiF.d_theta_L));
% Аппроксимация nλ = f(sin(θ)), неизвестный кристалл
unknown.fit2 = fit(sind(unknown.theta_L),unknown.n_L.*unknown.wl_L,{'x'});
% d (неизвестный кристалл), пм
unknown.d = coeffvalues(unknown.fit2)./2;
unknown.d_d = coeffvalues(unknown.fit2)./2;
%% Вывод результатов
printresult(wl_0,d_wl_0,'pm','λ','пм');
printresult(a,d_a,'pm','a','пм⋅кВ');
printresult(a_theor,0,'d','a_теор','пм⋅кВ');
printresult(LiF.wl_L,LiF.d_wl_L,'()','λ (L-серия)','пм');
printresult(2.*unknown.d,2.*unknown.d_d,'pm','2d','пм');
%% I(θ), разная сила тока эмиссии
figure('Name','I(θ), разная сила тока эмиссии','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:3
    plot(theta_J,I_J(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
legend('1 мА','0.8 мА','0.6 мА','Location','northeast');
pbaspect([2 1 1]);
figresize('width',17);
%% I(J)
figure('Name','I(J)','Units','centimeters');
hold on;

X = [0.6 1.0];
plot(X,polyval(polyfit(J,mean(I_J),1),X),'--','Color',[0.5 0.5 0.5]);
color = {'b','r','g'};
for i = 1:3
    plot(J(i),mean(I_J(:,i)),'o','Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itJ\rm, мА','FontName','Liberation Serif');
ylabel('\langle\itI\rm\rangle','FontName','Liberation Serif');
pbaspect([1 1 1]);
grid on;
figresize('width',9);
%% I(θ), характеристические линии, разная сила тока эмиссии
figure('Name','I(θ), характеристические линии, разная сила тока эмиссии','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:3
    plot(theta_ch_J,I_ch_J(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
legend('1 мА','0.8 мА','0.6 мА','Location','northeast');
pbaspect([4 3 1]);
figresize('width',17);
%% I(J), характеристические линии
figure('Name','I(J), характеристические линии','Units','centimeters');
hold on;

X = [0.6 1.0];
plot(X,polyval(polyfit(J,max(I_ch_J),1),X),'--','Color',[0.5 0.5 0.5]);
color = {'b','r','g'};
for i = 1:3
    plot(J(i),max(I_ch_J(:,i)),'o','Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itJ\rm, мА','FontName','Liberation Serif');
ylabel('\itI\rm_{max}','FontName','Liberation Serif');
pbaspect([1 1 1]);
grid on;
figresize('width',9);
%% I(θ), разное анодное напряжение и ток эмиссии
figure('Name','I(θ), разное анодное напряжение и ток эмиссии','Units','centimeters');
hold on;

plot(theta_J,I_J(:,3),'Color','g');
plot(theta_J,I_extra,'Color','k');

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
legend('35 кВ, 0.6 мА','27 кВ, 1 мА','Location','northeast');
pbaspect([2 1 1]);
figresize('width',17);
%% I(θ), разное анодное напряжение
figure('Name','I(θ), разное анодное напряжение','Units','centimeters');
hold on;

color = {'b','r','g','m'};
for i = 1:4
    plot(theta_U,I_U(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
legend('35 кВ','33 кВ','31 кВ','29 кВ','Location','northwest');
pbaspect([2 1 1]);
figresize('width',17);
%% λ_0(U^−1)
figure('Name','λ_0(U^−1)','Units','centimeters');
hold on;

X = [0.028 0.036];
plot(X,fit1(X),'--','Color',[0.5 0.5 0.5]);
color = {'b','r','g','m'};
for i = 1:4
    errorbar(U_inv(i),wl_0(i),d_wl_0(i),'o','Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm^{−1}, кВ^{−1}','FontName','Liberation Serif');
ylabel('\itλ\rm_0, пм','FontName','Liberation Serif');
pbaspect([1 1 1]);
grid on;
figresize('width',9);
%% I(θ), LiF
figure('Name','I(θ), LiF','Units','centimeters');
plot(LiF.theta,LiF.I,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
pbaspect([2 1 1]);
figresize('width',17);
%% I(θ), неизвестный кристалл
figure('Name','I(θ), неизвестный кристалл','Units','centimeters');
plot(unknown.theta,unknown.I,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('θ, °','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
pbaspect([2 1 1]);
figresize('width',17);
%% nλ = f(sin(θ)), неизвестный кристалл
figure('Name','nλ = f(sin(θ)), неизвестный кристалл','Units','centimeters');
hold on;

X = [0.4 0.8];
plot(X,unknown.fit2(X),'Color','g');
plot(sind(unknown.theta_L),unknown.n_L.*unknown.wl_L,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('sin(\itθ\rm)','FontName','Liberation Serif');
ylabel('\itnλ\rm, пм','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('height',9);