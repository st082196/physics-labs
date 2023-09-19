load('data_1');

% N
[N, d_N] = samplemean(sample_N);
% Аппроксимация линейного участка
fit_U_N = fit(U(fit_points),N(fit_points),{'x','1'},'Weights',1./d_N(fit_points).^2);
% N_ф
[N_noise, d_N_noise] = samplemean(sample_N_noise);
% M_i
[M_i, d_M_i] = samplemean(sample_M_i);
% m_i
m_i = M_i - N_noise;
d_m_i = sqrt(d_M_i.^2 + d_N_noise.^2);
% M_Σ
[M_sum, d_M_sum] = samplemean(sample_M_sum);
% m_Σ
m_sum = M_sum - N_noise;
d_m_sum = sqrt(d_M_sum.^2 + d_N_noise.^2);
% ∑m_i
sum_m_i = sources*m_i;
d_sum_m_i = sqrt(sources*d_m_i.^2);
% τ, мс
sum_m_i_sqr = sources*m_i.^2;
m_j = [repmat(m_i.',numel(m_sum),1),m_sum];
d_m_j = [repmat(d_m_i.',numel(m_sum),1),d_m_sum];
tau = (sum_m_i - m_sum)./(m_sum.^2 - sum_m_i_sqr).*t.*1e3;
d_tau = abs(tau).*sqrt(sum(((1./(sum_m_i-m_sum)+(2.*m_j)./(m_sum.^2-sum_m_i_sqr)).*d_m_j).^2.*[sources,ones(size(sources,1),1)],2));
p = 1./d_tau.^2;
mean_tau = sum(p.*tau)./sum(p);
d_mean_tau = 1./sum(p);
%% Вывод результатов
printresult(N,d_N,'pm','N');
printresult(N_noise,d_N_noise,'pm','N_ф');
printresult(M_i,d_M_i,'pm','M_i');
printresult(m_i,d_m_i,'pm','m_i');
printresult(M_sum,d_M_sum,'pm','M_Σ');
printresult(m_sum,d_m_sum,'pm','m_Σ');
printresult(sum_m_i,d_sum_m_i,'pm','∑m_i');
printresult(tau,d_tau,'pm','τ','мс');
printresult(mean_tau,d_mean_tau,'pm','mean_τ','мс');
%% N(U)
figure('Name','N(U)','Units','centimeters');
hold on;

X = U(fit_points);
plot(X,fit_U_N(X),'Color','g');
errorbar(U,N,d_N,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm, В','FontName','Liberation Serif');
ylabel('\itN','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);