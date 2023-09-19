load('data_2_altfit.mat');

[I, fit_N_I, gof_N_I, fit1, N_peaks, d_N, FWHM, R, S, S_0, F, E_graph, d_E_graph] = deal(cell(size(name)));
for i = 1:numel(name)
    I{i} = mean(sample_I{i},2);
    % Аппроксимация I(N)
    n_peaks = numel(N_peaks_approx{i});
    expression = 'a1*exp(-((x-b1)/c1)^2)';
    for j = 2:n_peaks
        expression = append(expression,sprintf(' + a%u*exp(-((x-b%u)/c%u)^2)',j,j,j));
    end
    expression = append(expression,' + d*x + e');
    [fit_N_I{i}, gof_N_I{i}] = fit( ...
        N_fit{i},I{i}(N_fit{i}),expression, ...
        'Robust','LAR', ...
        'StartPoint',[1000.*ones(n_peaks,1);N_peaks_approx{i};10.*ones(n_peaks,1);-1;0] ...
        );
    values = coeffvalues(fit_N_I{i});
    [fit1{i}.a, fit1{i}.b, fit1{i}.c] = matsplit(reshape(values(1:n_peaks*3),n_peaks,3),1);
    [fit1{i}.d, fit1{i}.e] = matsplit(values(end-1:end));
    variances = coeffvariances(fit_N_I{i});
    [fit1{i}.d_a, fit1{i}.d_b, fit1{i}.d_c] = matsplit(reshape(variances(1:n_peaks*3),n_peaks,3),1);
    [fit1{i}.d_d, fit1{i}.d_e] = matsplit(variances(end-1:end)); 
    [N_peaks{i}, order] = sort(fit1{i}.b);
    d_N{i} = fit1{i}.d_b;
    d_N{i} = d_N{i}(order);
    FWHM{i} = 2.*sqrt(log(2)).*fit1{i}.c;
    R{i} = FWHM{i}./N_peaks{i};
    S{i} = sqrt(pi).*fit1{i}.a.*fit1{i}.c;
    S_0{i} = trapz(I{i});
    F{i} = S{i}./S_0{i};
end
% Аппроксимация E(N)
N_cal = vertcat(N_peaks{nuclides_cal});
E_cal = vertcat(E_table{nuclides_cal});
weights = 1./vertcat(d_N{nuclides_cal}).^2;
fit_N_E = fit(N_cal,E_cal,{'x','1'},'Weights',weights);
[fit2.a, fit2.b] = matsplit(coeffvalues(fit_N_E));
[fit2.d_a, fit2.d_b] = matsplit(coeffvariances(fit_N_E));
% E_граф, кэВ
for i = 1:numel(name)
    E_graph{i} = fit_N_E(N_peaks{i});
    d_E_graph{i} = sqrt((fit2.a.*d_N{i}).^2 + ((N_peaks{i}-mean(N_cal)).*fit2.d_a).^2);
end
%% Вывод результатов
printresult(fit2.a,fit2.d_a,'pm','a','кэВ^−1');
printresult(fit2.b,fit2.d_b,'pm','b');
for i = 1:numel(name)
    fprintf('\n%s\n',name{i});
    printresult(N_peaks{i},d_N{i},'()','N');
    printresult(E_graph{i},d_E_graph{i},'()','E_граф','кэВ');
    printresult(R{i}.*100,1,'d','R','%');
    printresult(F{i}.*100,1,'d','F','%');
end
%% Вывод коэффициентов аппроксимации
i = 1;
printresult(fit1{i}.a,fit1{i}.d_a,'pm','a');
printresult(fit1{i}.b,fit1{i}.d_b,'pm','b');
printresult(fit1{i}.c,fit1{i}.d_c,'pm','c');
printresult(fit1{i}.d,fit1{i}.d_d,'pm','d');
printresult(fit1{i}.e,fit1{i}.d_e,'pm','e');
printresult(gof_N_I{i}.rmse,4,'s','RMSE');
%% I(N)
i = 1;

figure('Name','I(N)','Units','centimeters');
hold on;

plot(I{i},'Color','b');

set(gca,'box','on','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(name{i},'FontName','Liberation Serif');
xlabel('N','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
xlim([0 500]);
pbaspect([2 1 1]);
figresize('width',17);
%% I(N), аппроксимация
i = 1;

figure('Name','I(N)','Units','centimeters');
hold on;

x_fit = N_fit{i}(1):0.25:N_fit{i}(end);
for j=1:numel(N_peaks{i})
    plot(x_fit,fit1{i}.a(j).*exp(-((x_fit-fit1{i}.b(j))./fit1{i}.c(j)).^2)+fit1{i}.d.*x_fit+fit1{i}.e,'--','Color','g');
end
plot(x_fit,fit_N_I{i}(x_fit),'Color','g');
plot(N_fit{i},I{i}(N_fit{i}),'.','Color','b');
plot(remove(1:500,N_fit{i}),I{i}(remove(1:500,N_fit{i})),'x','Color','r','MarkerSize',2);

set(gca,'box','on','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title(name{i},'FontName','Liberation Serif');
xlabel('N','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
xlim([200 400]);
pbaspect([2 1 1]);
figresize('width',17);
%% E(N)
figure('Name','E(N)','Units','centimeters');
hold on;

x_fit = [0 350];
plot(x_fit,fit_N_E(x_fit),'Color','g');
plot(N_cal,E_cal,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('N','FontName','Liberation Serif','FontAngle','italic');
ylabel('E, кэВ','FontName','Liberation Serif','FontAngle','italic');
ylim([0 1600]);
grid on;
figresize('width',17);