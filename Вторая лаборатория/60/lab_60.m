load('data_60');

[U_at_extremum, d_U_at_extremum, I_at_extremum, d_I_at_extremum] = deal(cell(size(I)));
[delta_U, d_delta_U] = deal(NaN(max(cell2mat(section.'))-2,numel(section)));
for j = 1:numel(I)
    % Вычисление средних значений I(U) в экстремумах
    [U_at_extremum{j}, d_U_at_extremum{j}, I_at_extremum{j}, d_I_at_extremum{j}] = deal(zeros(max(section{j}),1));
    for i = 1:max(section{j})
        [U_at_extremum{j}(i), d_U_at_extremum{j}(i)] = samplemean(U{j}(section{j}==i),d_U_syst);
        [I_at_extremum{j}(i), d_I_at_extremum{j}(i)] = samplemean(I{j}(section{j}==i),d_I_syst);
    end
    % ∆U, В
    delta_U(1:numel(U_at_extremum{j})-2,j) = U_at_extremum{j}(3:end) - U_at_extremum{j}(1:end-2);
    d_delta_U(1:numel(U_at_extremum{j})-2,j) = sqrt(d_U_at_extremum{j}(1:end-2).^2 + d_U_at_extremum{j}(3:end).^2);
end
% ⟨∆U⟩, В
mean_delta_U = mean(delta_U,2,'omitnan');
d_mean_delta_U = sqrt(sum((d_delta_U).^2,2,'omitnan'))./sum(~isnan(d_delta_U),2);
%% Вывод результатов
for j = 1:numel(I)
    fprintf('\n\nU_к = %.1f В, I_н = %.2f А\n',U_c(j),I_f(j));
    printresult(U_at_extremum{j},d_U_at_extremum{j},'pm','U','В');
    printresult(I_at_extremum{j},d_I_at_extremum{j},'pm','I','нА');
    printresult(delta_U(~isnan(delta_U(:,j)),j),d_delta_U(~isnan(delta_U(:,j)),j),'pm','∆U','В');
end
printresult(mean_delta_U,d_mean_delta_U,'pm','⟨∆U⟩','В');
%% I(U)
figure('Name','I(U)','Units','centimeters');
hold on;

[x, y, X, Y] = deal(cell(size(I)));
marker = {'o','s','^'};
color = {'b','r','g'};
subset = zeros(size(color));
for j = 1:numel(color)
    x{j} = [U{j}(section{j}==0); U_at_extremum{j}];
    y{j} = [I{j}(section{j}==0); I_at_extremum{j}];

    X{j} = min(x{j}):0.1:max(x{j});
    Y{j} = pchip(x{j},y{j},X{j});

    plot(X{j},Y{j},'Color',color{j});
    subset(j) = plot(U{j}(section{j}==0),I{j}(section{j}==0),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3, ...
        'DisplayName',sprintf('\\itU\\rm_к = %.1f В,   \\itI\\rm_н = %.2f А',U_c(j),I_f(j)));
    plot(U_at_extremum{j},I_at_extremum{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',5);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm, В','FontName','Liberation Serif');
ylabel('\itI\rm, нА','FontName','Liberation Serif');
legend(subset,'Location','northwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);