load('data_2'); %#ok<*SAGROW>

[I, gof_N_I, N_peaks, d_N_peaks, FWHM, R, S, S_0, F, E.graph, d_E.graph] = deal(cell(size(name)));
for i = 1:numel(name)
    I{i} = mean(sample_I{i},2);
    % Аппроксимация I(N)
    num_peaks = numel(N_peaks_approx{i});
    [fit_N_I{i}.f, gof_N_I{i}.gof] = fit( ...
        fit_N_I{i}.points,I{i}(fit_N_I{i}.points),sprintf('gauss%u',num_peaks), ...
        'StartPoint',reshape([1000.*ones(num_peaks,1),N_peaks_approx{i},10.*ones(num_peaks,1)].',[],1) ...
        );
    [fit_N_I{i}.a, fit_N_I{i}.b, fit_N_I{i}.c] = matsplit(reshape(coeffvalues(fit_N_I{i}.f),3,num_peaks).',1);
    [fit_N_I{i}.d_a, fit_N_I{i}.d_b, fit_N_I{i}.d_c] = matsplit(reshape(coeffvariances(fit_N_I{i}.f),3,num_peaks).',1);
    [N_peaks{i}, order] = sort(fit_N_I{i}.b);
    d_N_peaks{i} = fit_N_I{i}.d_b(order);
    FWHM{i} = 2.*sqrt(log(2)).*fit_N_I{i}.c;
    R{i} = FWHM{i}./N_peaks{i};
    S{i} = sqrt(pi).*fit_N_I{i}.a.*fit_N_I{i}.c;
    S_0{i} = trapz(I{i});
    F{i} = S{i}./S_0{i};
end
% Аппроксимация E(N)
N_cal = vertcat(N_peaks{nuclides_cal});
E_cal = vertcat(E.ref{nuclides_cal});
fit_N_E.weights = 1./vertcat(d_N_peaks{nuclides_cal}).^2;
fit_N_E.f = fit(N_cal,E_cal,{'x','1'},'Weights',fit_N_E.weights);
[fit_N_E.a, fit_N_E.b] = matsplit(coeffvalues(fit_N_E.f)); % кэВ, кэВ
[fit_N_E.d_a, fit_N_E.d_b] = matsplit(coeffvariances(fit_N_E.f));
% E_граф, кэВ
for i = 1:numel(name)
    E.graph{i} = fit_N_E.f(N_peaks{i});
    d_E.graph{i} = sqrt((fit_N_E.a.*d_N_peaks{i}).^2 + ((N_peaks{i}-mean(N_cal)).*fit_N_E.d_a).^2);
end
%% Вывод результатов
printresult(fit_N_E.a,fit_N_E.d_a,'pm','a','кэВ');
printresult(fit_N_E.b,fit_N_E.d_b,'pm','b','кэВ');
printresult(vertcat(N_peaks{:}),vertcat(d_N_peaks{:}),'()','N');
printresult(vertcat(E.graph{:}),vertcat(d_E.graph{:}),'()','E_граф','кэВ');
printresult(vertcat(R{1:end}).*100,1,'d','R','%');
printresult(vertcat(F{1:end}).*100,1,'d','F','%');
%% I(N)
i = 1;

figure('Name','I(N)','Units','centimeters');
plot(I{i},'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
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

X = (1:0.5:500).';
Y = zeros(numel(X),numel(N_peaks{i}));
for j = 1:numel(N_peaks{i})
    Y(:,j) = fit_N_I{i}.a(j).*exp(-((X-fit_N_I{i}.b(j))./fit_N_I{i}.c(j)).^2);
    plot(X,Y(:,j),'--','Color','g');
end
plot(X,fit_N_I{i}.f(X),'Color','g');
plot(fit_N_I{i}.points,I{i}(fit_N_I{i}.points),'.','Color','b');

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
title(name{i},'FontName','Liberation Serif');
xlabel('N','FontName','Liberation Serif','FontAngle','italic');
ylabel('I','FontName','Liberation Serif','FontAngle','italic');
xlim([0 500]);
pbaspect([2 1 1]);
figresize('width',17);
%% E(N)
figure('Name','E(N)','Units','centimeters');
hold on;

X = [0 350];
plot(X,fit_N_E.f(X),'Color','g');
plot(N_cal,E_cal,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itN','FontName','Liberation Serif');
ylabel('\itE\rm, кэВ','FontName','Liberation Serif');
ylim([0 1600]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);