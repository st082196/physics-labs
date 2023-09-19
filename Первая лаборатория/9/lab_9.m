load('data_9');
%% 1.
unpackstruct(part1);

% Аппроксимация S(t)
fit1 = fit(t,S,{'x','1'});
[c.exp, S_0] = matsplit(coeffvalues(fit1).*[1e1 1]); % м/с, см
[d_c.exp, d_S_0] = matsplit(coeffvariances(fit1).*[1e1 1]);
% c_теор, м/с
c.theor = sqrt(gamma.theor.*R.*(theta+273.15)./(M.*1e-3));
d_c.theor = abs(c.theor./sqrt(2).*d_theta./(theta+273.15));
%% Вывод результатов
printresult(c.exp,d_c.exp,'pm','c_эксп','м/с');
printresult(S_0,d_S_0,'pm','S_0','см');
printresult(c.theor,d_c.theor,'pm','c_теор','м/с');
%% S(t)
figure('Name','S(t)','Units','centimeters');
hold on;

X = [0.4 1.1];
plot(X,fit1(X),'Color','g');
plot(t,S,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itt\rm, мс','FontName','Liberation Serif');
ylabel('\itS\rm, см','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% c, м/с
c.exp = (S - S_0)./t.*1e1;
func_c = @(theta) sqrt(gamma.theor.*R.*(theta+273.15)./(M.*1e-3));
c.theor = func_c(theta);
% ∆c, %
delta_c = (c.exp - c.theor)./c.theor.*100;
% Аппроксимация c^2(θ)
fit2 = fit(theta,c.exp.^2,{'x','1'});
[a, b] = matsplit(coeffvalues(fit2)); % м^2⋅с^−2⋅°C^−1, м^2⋅с^−2
[d_a, d_b] = matsplit(coeffvariances(fit2));
% γ_эксп
gamma.exp = M./R.*a.*1e-3;
d_gamma.exp = abs((M.*1e-3)./R.*d_a);
%% Вывод результатов
printresult(c.exp,4,'s','c_эксп','м/с');
printresult(c.theor,4,'s','c_теор','м/с');
printresult(delta_c,3,'s','∆c','%');
printresult(a,d_a,'pm','a','м^2⋅с^−2⋅°C^−1','auto');
printresult(b,d_b,'pm','b','м^2⋅с^−2','auto');
printresult(gamma.exp,d_gamma.exp,'pm','γ_эксп');
%% c(θ)
figure('Name','c(θ)','Units','centimeters');
hold on;

X = 20:1:60;
plot(X,func_c(X),'Color','g');
plot(theta,c.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itθ\rm, °C','FontName','Liberation Serif');
ylabel('\itc\rm, м/с','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% c^2(θ)
figure('Name','c^2(θ)','Units','centimeters');
hold on;

X = [20 60];
plot(X,func_c(X).^2,'Color','g');
plot(theta,c.exp.^2,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(X,fit2(X),'Color','b');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itθ\rm, °C','FontName','Liberation Serif');
ylabel('\itc\rm^2, м^2/с^2','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);