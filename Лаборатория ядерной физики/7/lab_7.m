load('data_7');

% Аппроксимация E(N)
N_cal = N_peaks{1}(peaks_cal);
E_cal = E_peaks.ref{1}(peaks_cal);
weights = 1./d_N_peaks{1}(peaks_cal).^2;
fit_N_E = fit(N_cal,E_cal,{'x','1'},'Weights',weights);
[a, b] = matsplit(coeffvalues(fit_N_E)); % кэВ, кэВ
if numel(peaks_cal) > 2
    [d_a, d_b] = matsplit(coeffvariances(fit_N_E));
else
    [d_a, d_b] = deal(0);
end

N = (1:numel(I{1})).';
E = fit_N_E(N);
[E_peaks.exp, d_E_peaks.exp, delta_N_peaks, delta_E_peaks, R] = deal(cell(4,1));
for i = 1:numel(name)
    % E, кэВ
    E_peaks.exp{i} = fit_N_E(N_peaks{i});
    d_E_peaks.exp{i} = sqrt((a.*d_N_peaks{i}).^2 + ((N_peaks{i}-mean(N_cal)).*d_a).^2);
    % ∆N
    for peak = 1:numel(N_peaks{i})
        peak_N = peak_bounds{i}(peak,1):peak_bounds{i}(peak,2);
        delta_N_peaks{i}(peak,1) = peakwidth(peak_N,I{i}(peak_N));
    end
    % ∆E, кэВ
    delta_E_peaks{i} = delta_N_peaks{i}.*a;
    % R, %
    R{i} = delta_E_peaks{i}./E_peaks.exp{i}.*100;
end
I_1 = trapz(I{1}(peak_bounds{1}(1,1):peak_bounds{1}(1,2))) + trapz(I{1}(peak_bounds{1}(2,1):peak_bounds{1}(2,2)));
I_5 = trapz(I{1}(peak_bounds{1}(3,1):peak_bounds{1}(3,2)));
% λ, с^−1
lambda = log(2)./T;
% t, дн
c = zeros(numel(lambda),1);
for i = 1:numel(c)
    c(i) = prod(lambda)./prod(remove(lambda,i) - lambda(i));
end
f = @(t) sum(c.*exp(-lambda.*t))./lambda(1).*exp(-lambda(1).*t) - I_5./I_1;
t = fzero(f,1e7)./86400;
%% Вывод результатов
printresult(a,d_a,'pm','a','кэВ');
printresult(b,d_b,'pm','b','кэВ');
printresult(vertcat(N_peaks{:}),vertcat(d_N_peaks{:}),'pm','N');
printresult(vertcat(delta_N_peaks{:}),0,'d','∆N');
printresult(vertcat(E_peaks.exp{:}),vertcat(d_E_peaks.exp{:}),'pm','E','кэВ');
printresult(vertcat(delta_E_peaks{:}),0,'d','∆E','кэВ',0);
printresult(vertcat(R{:}),2,'d','R','%');
printresult(I_1,0,'d','I_1');
printresult(I_5,0,'d','I_5');
printresult(t,0,'d','t','дн');
%% I(N)
i = 1;

figure('Name','I(N)','Units','centimeters');
plot(I{i},'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(name{i},'FontName','Liberation Serif');
xlabel('N','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
xlim([0 N(end)]);
pbaspect([2 1 1]);
figresize('width',17);
%% E(N)
figure('Name','E(N)','Units','centimeters');
hold on;

X = [2000 3600];
plot(X,fit_N_E(X),'Color','g');
plot(N_cal,E_cal,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
% plot(vertcat(N_peaks{2:end}),vertcat(E_peaks.ref{2:end}),'o','Color','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itN','FontName','Liberation Serif');
ylabel('\itE\rm, кэВ','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% I(E)
i = 1;

figure('Name','I(E)','Units','centimeters');
plot(E,I{i},'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(name{i},'FontName','Liberation Serif');
xlabel('\itE\rm, кэВ','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
xlim([0 E(end)]);
pbaspect([2 1 1]);
figresize('width',17);
%% f(t)
figure('Name','f(t)');
x = logspace(0,8,1e6);
semilogx(x,f(x));
yline(0);