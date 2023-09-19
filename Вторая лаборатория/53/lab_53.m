load('data_53');
%% 1.
unpackstruct(part1);

% Аппроксимация I(U_с)
fits = struct('f',cell(size(I)));
best_fit = struct('adjrsquare',cell(size(I)));
[R, d_R] = deal(zeros(size(I)));
for j = 1:numel(I)
    [fits(j).f, fits(j).gof] = deal(cell(numel(I{j})-min_num_fit_points+1,1));
    for n = min_num_fit_points:numel(I{j})
        i = n - min_num_fit_points + 1;
        [fits(j).f{i}, fits(j).gof{i}] = fit(U_d{j}(1:n),I{j}(1:n),{'x'});
    end
    fits(j).gof = [fits(j).gof{:}].';
    [best_fit(j).adjrsquare, best_fit(j).i] = max([fits(j).gof.adjrsquare]);
    best_fit(j).n = best_fit(j).i + min_num_fit_points - 1;
    best_fit(j).x = U_d{j}(1:best_fit(j).n);
    best_fit(j).y = I{j}(1:best_fit(j).n);
    best_fit(j).f = fits(j).f{best_fit(j).i};
    best_fit(j).a = coeffvalues(best_fit(j).f); % мА/В
    best_fit(j).d_a = coeffvariances(best_fit(j).f);
    % R, Ом
    R(j) = 1./best_fit(j).a.*1e3;
    d_R(j) = abs(best_fit(j).d_a./best_fit(j).a.^2).*1e3;
end
%% Вывод результатов
printresult([best_fit.a],[best_fit.d_a],'pm','a','мА/В');
printresult(R,d_R,'pm','R','Ом');
%% I(U_с)
figure('Name','I(U_с)','Units','centimeters');
hold on;

[X, Y] = deal(cell(size(I)));
marker = {'o','s','^','v','p'};
color = {'b','r','g','m','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    X{j} = (min(U_d{j}):0.01:max(U_d{j})).';
    Y{j} = pchip(U_d{j},I{j},X{j});
    plot(X{j},Y{j},'Color',color{j});
    subset(j) = plot(U_d{j},I{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3, ...
        'DisplayName',sprintf('\\itU\\rm_з = %+.1f В',U_g(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_с, В','FontName','Liberation Serif');
ylabel('\itI\rm, мА','FontName','Liberation Serif');
legend(subset,'Location','southeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% I(U_с), линейный участок
figure('Name','I(U_с), линейный участок','Units','centimeters');
hold on;

X = [0 0.7];
marker = {'o','s','^','v','p'};
color = {'b','r','g','m','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(X,best_fit(j).f(X),'Color',color{j});
    subset(j) = plot(best_fit(j).x,best_fit(j).y,marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3, ...
        'DisplayName',sprintf('\\itU\\rm_з = %+.1f В',U_g(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('Линейный участок');
xlabel('\itU\rm_с, В','FontName','Liberation Serif');
ylabel('\itI\rm, мА','FontName','Liberation Serif');
legend(subset,'Location','southeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% Аппроксимация I(U_з)
fits = struct('f',cell(size(I)));
best_fit = struct('adjrsquare',cell(size(I)));
[U_g_cutoff, d_U_g_cutoff] = deal(zeros(size(I)));
for j = 1:numel(I)
    [fits(j).f, fits(j).gof] = deal(cell(numel(I{j})-min_num_fit_points+1,1));
    for n = min_num_fit_points:numel(I{j})
        i = n - min_num_fit_points + 1;
        [fits(j).f{i}, fits(j).gof{i}] = fit(U_g{j}(1:n),I{j}(1:n),{'x','1'});
    end
    fits(j).gof = [fits(j).gof{:}].';
    [best_fit(j).adjrsquare, best_fit(j).i] = max([fits(j).gof.adjrsquare]);
    best_fit(j).n = best_fit(j).i + min_num_fit_points - 1;
    best_fit(j).x = U_g{j}(1:best_fit(j).n);
    best_fit(j).y = I{j}(1:best_fit(j).n);
    best_fit(j).f = fits(j).f{best_fit(j).i};
    [best_fit(j).a, best_fit(j).b] = matsplit(coeffvalues(best_fit(j).f)); % мА/В, мА
    [best_fit(j).d_a, best_fit(j).d_b] = matsplit(coeffvariances(best_fit(j).f));
    % U_зо, В
    U_g_cutoff(j) = -best_fit(j).b./best_fit(j).a;
    d_U_g_cutoff(j) = abs(mean(best_fit(j).y).*best_fit(j).d_a./best_fit(j).a.^2);
end
%% Вывод результатов
printresult([best_fit.a],[best_fit.d_a],'pm','a','мА/В');
printresult([best_fit.b],[best_fit.d_b],'pm','b','мА');
printresult(U_g_cutoff,d_U_g_cutoff,'pm','U_зо','В');
%% I(U_з)
figure('Name','I(U_з)','Units','centimeters');
hold on;

[X, Y] = deal(cell(size(I)));
marker = {'o','s','^','v','p'};
color = {'b','r','g','m','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    X{j} = (min(U_g{j}):0.01:max(U_g{j})).';
    Y{j} = pchip(U_g{j},I{j},X{j});
    plot(X{j},Y{j},'Color',color{j});
    subset(j) = plot(U_g{j},I{j},marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3, ...
        'DisplayName',sprintf('\\itU\\rm_с = %.1f В',U_d(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_з, В','FontName','Liberation Serif');
ylabel('\itI\rm, мА','FontName','Liberation Serif');
legend(flip(subset),'Location','northwest','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% I(U_з), линейный участок
figure('Name','I(U_з), линейный участок','Units','centimeters');
hold on;

marker = {'o','s','^','v','p'};
color = {'b','r','g','m','k'};
subset = zeros(size(color));
for j = 1:numel(color)
    X = [U_g_cutoff(j) 0];
    plot(X,fits(j).f{best_fit(j).i}(X),'Color',color{j});
    subset(j) = plot(best_fit(j).x,best_fit(j).y,marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3, ...
        'DisplayName',sprintf('\\itU\\rm_с = %.1f В',U_d(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('Линейный участок');
xlabel('\itU\rm_з, В','FontName','Liberation Serif');
ylabel('\itI\rm, мА','FontName','Liberation Serif');
legend(flip(subset),'Location','northwest','FontName','Liberation Serif');
xlim([-2 0]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);