load('data_114');

% x, см
[x_1, d_x_1] = samplemean(sample_x_1,d_x_syst);
[x_2, d_x_2] = samplemean(sample_x_2,d_x_syst);
% R, мм
[R_1, d_R_1] = samplemean(sample_R_1,d_R_syst);
[R_2, d_R_2] = samplemean(sample_R_2,d_R_syst);
% d, см
d = abs(x_2 - x_1);
d_d = sqrt(d_x_1.^2 + d_x_2.^2);
% L, см
L = d.*(sqrt(R_1) + sqrt(R_2))./(sqrt(R_1) - sqrt(R_2));
d_L = 1./(sqrt(R_1.*R_2).*(sqrt(R_1) - sqrt(R_2)).^2) .* ...
    sqrt((sqrt(R_1.*R_2).*(R_1 - R_2).*d_d).^2 + (d.*R_2.*d_R_1).^2 + (d.*R_1.*d_R_2).^2);
% 2t, мм
two_t = sqrt(R_1.*R_2);
d_two_t = 1./(2.*sqrt(R_1.*R_2)).*sqrt((R_2.*d_R_1).^2 + (R_1.*d_R_2).^2);
% Аппроксимация y(m)
fit1 = cell(size(y));
fit2 = struct('x',repmat({L./two_t.*1e1},1,2));
[sigma, y_0, d_sigma, d_y_0, fit2(1).y, fit2(2).y] = deal(zeros(size(y)));
for i = 1:numel(y)
    fit1{i} = fit(m,y{i},{'x','1'});
    [sigma(i), y_0(i)] = matsplit(coeffvalues(fit1{i})); % мм, мм
    [d_sigma(i), d_y_0(i)] = matsplit(coeffvariances(fit1{i}));
    fit2(1).y(i) = sigma(i);
    fit2(2).y(i) = y{i}(m == 0.5) - y{i}(m == -0.5);
end
% Аппроксимация σ(L/2t)
for i = 1:numel(fit2)
    fit2(i).f = fit(fit2(i).x,fit2(i).y,{'x'});
    % λ, нм
    fit2(i).wl = coeffvalues(fit2(i).f).*1e6;
    fit2(i).d_wl = coeffvariances(fit2(i).f).*1e6;
end
%% Вывод результатов
printresult(x_1,d_x_1,'pm','x_1','см');
printresult(x_2,d_x_2,'pm','x_2','см');
printresult(R_1,d_R_1,'pm','R_1','мм');
printresult(R_2,d_R_2,'pm','R_2','мм');
printresult(d,d_d,'pm','d','см');
printresult(L,d_L,'pm','L','см');
printresult(two_t,d_two_t,'pm','2t','мм');
printresult(sigma,d_sigma,'pm','σ','мм');
printresult([fit2.wl],[fit2.d_wl],'pm','λ','нм');
%% y(m)
figure('Name','y(m)','Units','centimeters');
hold on;

x_fit = m([1 end]);
marker = {'o','s','^'};
color = {'b','r','g'};
subset = zeros(size(color));
for i = 1:numel(color)
    plot(x_fit,fit1{i}(x_fit)-y_0(i),'Color',color{i});
    subset(i) = plot(m,y{i}-y_0(i),marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itm','FontName','Liberation Serif');
ylabel('\ity\rm − \ity\rm_0, мм','FontName','Liberation Serif');
legend(subset,'\itx\rm_0 = 55 см','\itx\rm_0 = 60 см','\itx\rm_0 = 65 см','Location','southeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% σ(L/2t)
i = 1;

figure('Name','σ(L/2t)','Units','centimeters');
hold on;

x_fit = [540 680];
plot(x_fit,fit2(i).f(x_fit),'Color','g');
plot(fit2(i).x,fit2(i).y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itL\rm/2\itt','FontName','Liberation Serif');
switch i
    case 1
        ylabel('\itσ\rm, мм','FontName','Liberation Serif');
    case 2
        ylabel('\itσ\rm = \ity\rm_{0.5} − \ity\rm_{−0.5}, мм','FontName','Liberation Serif');
end
pbaspect([4 3 1]);
grid on;
figresize('width',17);