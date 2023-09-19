load('data_5');

% ∆_N
d_N = sqrt(abs(N)).*d_N_1;
% ∆n
delta_n = N.*wl./L.*1e-9;
d_delta_n = abs(delta_n).*sqrt((d_N./N).^2 + (d_wl./wl).^2 + (d_L./L).^2);
d_delta_n(delta_n==0) = 0;
% ∆p, кПа
delta_p = rho.*g.*(h_a-h_c).*1e-6;
d_delta_p = abs(delta_p).*sqrt((d_rho./rho).^2 + (d_g./g).^2 + 2.*(d_h./(h_a-h_c)).^2);
d_delta_p(delta_p==0) = sqrt(2).*abs(rho.*g.*d_h).*1e-6;
% Аппроксимация ∆n(∆p)
fit1 = fit(reshape(delta_p,[],1),repmat(delta_n,size(delta_p,2),1),{'x'});
% β, кПа^−1
beta = coeffvalues(fit1);
d_beta = coeffvariances(fit1);
% n
n = 1 + beta.*p;
d_n = sqrt((p.*0.02.*beta).^2 + (beta.*d_p).^2);
%% Вывод результатов
printresult(N(N~=0),d_N(N~=0),'pm','N');
printresult(delta_n(N~=0),d_delta_n(N~=0),'pm','∆n',[],'auto');
printresult(delta_p(:,1),d_delta_p(:,1),'pm','Серия №1, ∆p','кПа');
printresult(delta_p(:,2),d_delta_p(:,2),'pm','Серия №2, ∆p','кПа');
printresult(beta,d_beta,'pm','β','кПа^−1','auto');
printresult(n,d_n,'()','n');
%% ∆n(∆p)
figure('Name','∆n(∆p)','Units','centimeters');
hold on;

X = [-2 4];
plot(X,fit1(X),'Color','g');
plot(delta_p(:,1),delta_n,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(delta_p(:,2),delta_n,'s','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\it∆p\rm, кПа','FontName','Liberation Serif');
ylabel('\it∆n','FontName','Liberation Serif');
legend('','Серия №1','Серия №2','Location','southeast');
ylim([-5 10].*1e-6);
pbaspect([4 3 1]);
grid on;
figresize('height',10);