load('data_8');

% N_ф, имп/с
[N_noise, d_N_noise] = samplemean(sample_N_noise);
% N′, имп/с
[N_prime, d_N_prime] = samplemean(sample_N_prime);
% N, имп/с
N = N_prime - N_noise;
d_N = sqrt(d_N_prime.^2 + d_N_noise.^2);
% Аппроксимация N(x)
[fit_x_N, gof_x_N] = fit(x,N,'a/(x-x0)^2','Weights',1./d_N.^2,'Robust','LAR','StartPoint',[1e4 0]);
[a, x_0] = matsplit(coeffvalues(fit_x_N)); % имп⋅c^−1⋅см^2, см
[d_a, d_x_0] = matsplit(coeffvariances(fit_x_N));
% A, Бк
A = 16.*a./(D.^2.*eps).*1e2;
d_A = abs(A.*d_a./a);
% d, мм
[d, d_d] = samplemean(sample_d,d_d_syst);
% M_ф, имп/с
[M_noise, d_M_noise] = samplemean(sample_M_noise);
% M′, имп/с
[M_prime, d_M_prime] = samplemean(sample_M_prime);
% M, имп/с
M = M_prime - M_noise;
d_M = sqrt(d_M_prime.^2 + d_M_noise.^2);
% N/M
N_div_M = N(1)./M;
d_N_div_M = abs(N_div_M).*sqrt((d_N(1)./N(1)).^2 + (d_M./M).^2);
% μ, см^−1
mu = log(N(1)./M)./d.*10;
d_mu = abs(1./d).*sqrt((d_N(1)./N(1)).^2 + (d_M./M).^2 + (log(N(1)./M)./d.*d_d).^2).*10;
% P_ф, мкЗв/ч
[P_noise, d_P_noise] = samplemean(sample_P_noise,d_P_syst);
% P′, мкЗв/ч
[P_prime, d_P_prime] = samplemean(sample_P_prime,d_P_syst);
% P, мкЗв/ч
P = P_prime - P_noise;
d_P = sqrt(d_P_prime.^2 + d_P_noise.^2);
% P_расч, мкЗв/ч
P_calc = K_gamma.*A./R.^2.*1e-3;
d_P_calc = abs(P_calc).*sqrt((d_A./A).^2 + (2.*d_R./R).^2);
% Эффективная доза, мкЗв
E = R.^2./R_E.^2.*P.*t;
% N′_1, имп/с
[N_prime_1, d_N_prime_1] = samplemean(sample_N_prime_1);
% N_1, имп/с
N_1 = N_prime_1 - N_noise_1;
d_N_1 = d_N_prime_1;
% N′_2, имп/с
[N_prime_2, d_N_prime_2] = samplemean(sample_N_prime_2);
% N_2, имп/с
N_2 = N_prime_2 - N_noise_2;
d_N_2 = d_N_prime_2;
% N_совп, имп/с
N_coinc = N_prime_coinc - N_noise_coinc;
% N_сл, имп/с
N_rand = 2.*tau.*N_1.*N_2.*1e-6;
% A (метод совпадений), Бк
A_coinc = N_1.*N_2./(N_coinc - N_rand);
d_A_coinc = abs(A_coinc).*sqrt((d_N_1./N_1).^2 + (d_N_2./N_2).^2);
%% Вывод результатов
printresult(N_noise,d_N_noise,'pm','N_ф','имп/с');
printresult(N_prime,d_N_prime,'pm','N′','имп/с');
printresult(N,d_N,'pm','N','имп/с');
printresult(a,d_a,'pm','a','имп⋅c^−1⋅см^2');
printresult(x_0,d_x_0,'pm','x_0','см');
printresult(A,d_A,'pm','A','Бк');
printresult(A.*Bq2Ci,d_A.*Bq2Ci,'pm','A','Ки','auto');
printresult(d,d_d,'pm','d','мм');
printresult(M_noise,d_M_noise,'pm','M_ф','имп/с');
printresult(M_prime,d_M_prime,'pm','M′','имп/с');
printresult(M,d_M,'pm','M','имп/с');
printresult(N_div_M,d_N_div_M,'pm','N/M');
printresult(mu,d_mu,'pm','μ','см^−1');
printresult(P_noise,d_P_noise,'pm','P_ф','мкЗв/ч');
printresult(P_prime,d_P_prime,'pm','P′','мкЗв/ч');
printresult(P,d_P,'pm','P','мкЗв/ч');
printresult(P_calc,d_P_calc,'pm','P_расч','мкЗв/ч');
printresult(E(1,:),2,'s','Эффективная доза за 6 ч','мкЗв');
printresult(E(2,:),2,'s','Эффективная доза за год','мкЗв');
printresult(N_prime_1,d_N_prime_1,'pm','N′_1','имп/с');
printresult(N_1,d_N_1,'pm','N_1','имп/с');
printresult(N_prime_2,d_N_prime_2,'pm','N′_2','имп/с');
printresult(N_2,d_N_2,'pm','N_2','имп/с');
printresult(N_coinc,3,'s','N_совп','имп/с');
printresult(N_rand,2,'s','N_сл','имп/с');
printresult(A_coinc,d_A_coinc,'pm','A','Бк');
printresult(A_coinc.*Bq2Ci,d_A_coinc.*Bq2Ci,'pm','A','Ки','auto');
%% N(x)
figure('Name','N(x)','Units','centimeters');
hold on;

X = 12:0.1:26;
plot(X,fit_x_N(X),'Color','g');
plot(x,N,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx\rm, см','FontName','Liberation Serif');
ylabel('\itN','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);