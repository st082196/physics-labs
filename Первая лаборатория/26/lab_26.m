load('data_26');
%% 1.
unpackstruct(part1);

% Аппроксимация U(I)
fit1 = cell(size(U));
labels = strings(size(fit1));
[R, d_R] = deal(zeros(size(fit1)));
for i = 1:numel(fit1)
    fit1{i} = fit(I{i},U{i},{'x'});
    labels(i) = sprintf('%s\\fontname{Liberation Serif}, \\itd\\rm = %.2f мм, \\itl\\rm = %.2f м',material(i),d(i),l(i));
    % R, Ом
    R(i) = coeffvalues(fit1{i});
    d_R(i) = coeffvariances(fit1{i});
end
% S, мм^2
S = pi.*d.^2./4;
d_S = abs(pi.*d./2.*d_d);
% ρ, Ом⋅мм^2/м
rho = R.*pi.*d.^2./(4.*l);
d_rho = abs(rho).*sqrt((d_R./R).^2 + (2.*d_d./d).^2 + (d_l./l).^2);
[mean_rho, d_mean_rho] = deal(zeros(size(unique(material))));
i = 0;
for material_i = unique(material)
    i = i + 1;
    points = material == material_i;
    mean_rho(i) = mean(rho(points));
    d_mean_rho(i) = norm(d_rho(points))./nnz(points);
end
% Аппроксимация U(I), d = const
fit2.f = fit([I{2};I{5}],[U{2}.*2;U{5}],{'x'});
fit2.R = coeffvalues(fit2.f); % Ом
fit2.d_R = coeffvariances(fit2.f);
% Аппроксимация R(S^−1), l = const
fit3.x = S(1:4).'.^-1;
fit3.y = R(1:4).';
fit3.f = fit(fit3.x,fit3.y,{'x'});
%% Вывод результатов
printresult(S,d_S,'pm','S','мм^2');
printresult(R,d_R,'pm','R','Ом');
printresult(rho,d_rho,'pm','ρ','Ом⋅мм^2/м');
printresult(mean_rho,d_mean_rho,'pm','ρ (ср.)','Ом⋅мм^2/м');
printresult(fit2.R,fit2.d_R,'pm','R (d = const)','Ом');
printresult(fit2.R./2,fit2.d_R./2,'pm','R/2 (d = const)','Ом');
%% U(I)
figure('Name','U(I)','Units','centimeters');
hold on;

x_fit = [0 2];
marker = {'o','s','^','v','d','p'};
color = {'b','r','g','m','c','k'};
subset = zeros(size(color));
for i = 1:numel(color)
    plot(x_fit,fit1{i}(x_fit),'Color',color{i});
    subset(i) = plot(I{i},U{i},marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itI\rm, А','FontName','Liberation Serif');
ylabel('\itU\rm, В','FontName','Liberation Serif');
legend(subset,labels,'Location','northwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% U(I), d = const
figure('Name','U(I), d = const','Units','centimeters');
hold on;

x_fit = [0 2];
plot(x_fit,fit2.f(x_fit),'Color','g');
plot(I{2},U{2}.*2,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
plot(I{5},U{5},'s','Color','r','MarkerFaceColor','r','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('\itd\rm\bf = const','FontName','Liberation Serif');
xlabel('\itI\rm, А','FontName','Liberation Serif');
ylabel('\itU\rm, В','FontName','Liberation Serif');
legend('',labels(2)+' (2\itU\rm)',labels(5),'Location','northwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% R(S^−1), l = const
figure('Name','R(S^−1), l = const','Units','centimeters');
hold on;

x_fit = [0 12];
plot(x_fit,fit3.f(x_fit),'Color','g');
plot(fit3.x,fit3.y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('\itl\rm\bf = const','FontName','Liberation Serif');
xlabel('\itS\rm^{−1}, мм^{−2}','FontName','Liberation Serif');
ylabel('\itR\rm, Ом','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% Аппроксимация U(I)
fit4 = fit(I,U,{'x','1'});
[r, U_0] = matsplit(coeffvalues(fit4).*[-1,1]); %
[d_r, d_U_0] = matsplit(coeffvariances(fit4));
% I_кз, А
I_sc = U_0./r;
d_I_sc = abs(I_sc).*sqrt((d_U_0./U_0).^2 + (d_r./r).^2);
% P_2, max, 
P_2_max = I_sc.*U_0;
d_P_2_max = sqrt((U_0.*d_I_sc).^2 + (I_sc.*d_U_0).^2);
% R, Ом
R = U./I;
% I^*, I_кз
I_r = I./I_sc;
% β, r
beta = R./r;
% P, Вт
P_1 = I.*U;
P_2 = I.*U_0;
% P^*, P_2, max
P_1_r = P_1./P_2_max;
P_2_r = P_2./P_2_max;
% η
eta = P_1./P_2;
%% Вывод результатов
printresult(r,d_r,'pm','r','Ом');
printresult(U_0,d_U_0,'pm','U_0','В');
printresult(I_sc,d_I_sc,'pm','I_кз','А');
printresult(P_2_max,d_P_2_max,'pm','P_2, max','Вт');
printresult(R,2,'s','R','Ом');
printresult(I_r,2,'s','I^*','I_кз');
printresult(beta,2,'s','β','r');
printresult(P_1_r,2,'s','P_1^*','P_2, max');
printresult(P_2_r,2,'s','P_2^*','P_2, max');
printresult(eta.*100,2,'s','η','%');
%% U(I)
figure('Name','U(I)','Units','centimeters');
hold on;

x_fit = [0 0.6];
plot(x_fit,fit4(x_fit),'Color','g');
plot(I,U,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itI\rm, А','FontName','Liberation Serif');
ylabel('\itU\rm, В','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% P_1^*(I^*), P_2^*(I^*), η(I^*)
figure('Name','P_1^*(I^*), P_2^*(I^*), η(I^*)','Units','centimeters');
hold on;

plot(I_r,P_1_r,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'DisplayName','\itP\rm_1^*, \itP\rm_{2, max}');
plot(I_r,P_2_r,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'DisplayName','\itP\rm_2^*, \itP\rm_{2, max}');
plot(I_r,eta,'^','Color','g','MarkerFaceColor','g','MarkerSize',3,'DisplayName','\itη');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itI\rm^*, \itI\rm_{кз}','FontName','Liberation Serif');
legend('Location','east','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% P_1^*(β), P_2^*(β), η(β)
figure('Name','P_1^*(β), P_2^*(β), η(β)','Units','centimeters');
hold on;

plot(beta,P_1_r,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'DisplayName','\itP\rm_1^*, \itP\rm_{2, max}');
plot(beta,P_2_r,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'DisplayName','\itP\rm_2^*, \itP\rm_{2, max}');
plot(beta,eta,'^','Color','g','MarkerFaceColor','g','MarkerSize',3,'DisplayName','\itη');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itβ\rm, \itr','FontName','Liberation Serif');
legend('Location','east','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);