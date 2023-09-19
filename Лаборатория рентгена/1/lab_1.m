load('data_1');
%% 1.
unpackstruct(part1);

T = I(2:end,:)./I(1,:);
% Аппроксимация ln(T) = f(d)
[x, y, fit1] = deal(cell(1,size(T,2)));
[mu, d_mu] = deal(zeros(size(fit1)));
for i = 1:numel(fit1)
    x{i} = d;
    y{i} = log(T(:,i));
    fit1{i} = fit(x{i},y{i},{'x'});
    % μ, см^−1
    mu(i) = -coeffvalues(fit1{i}).*1e1;
    d_mu(i) = coeffvariances(fit1{i}).*1e1;
end
% μ_m, см^2/г
mu_m = mu./rho;
d_mu_m = abs(d_mu./rho);
% σ_m, см^2/г
sigma_m = 8.*pi./3.*r_e.^2.*Z.*N_A./M.*1e4;
% τ_m, см^2/г
tau_m = mu_m - sigma_m;
d_tau_m = d_mu_m;
%% Вывод результатов
printresult(T(:,1),5,'d','T без фильтра');
printresult(T(:,2),5,'d','T с фильтром');
printresult(mu,d_mu,'pm','μ','см^−1');
printresult(mu_m,d_mu_m,'pm','μ_m','см^2/г');
printresult(sigma_m,2,'s','σ_m','см^2/г');
printresult(tau_m,d_tau_m,'pm','τ_m','см^2/г');
%% ln(T) = f(d)
figure('Name','ln(T) = f(d)','Units','centimeters');
hold on;

X = [0 3];
marker = {'o','s'};
color = {'b','r'};
subset = zeros(size(color));
for i = 1:numel(color)
    plot(X,fit1{i}(X),'Color',color{i});
    subset(i) = plot(x{i},y{i},marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itd\rm, мм','FontName','Liberation Serif');
ylabel('ln(\itT\rm)','FontName','Liberation Serif');
legend(subset,'без фильтра','с фильтром','Location','northeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% Приведение интенсивности к одинаковому эмиссионному току
I(1:3,:) = I(1:3,:).*50;
T = I(2:end,:)./I(1,:);
% μ, см^−1
mu = -log(T)./d.*10;
Z_K = 1./sqrt(wl_K.*1e-12.*R);
%% Вывод результатов
fprintf('\nБез фильтра\n');
printresult(T(:,1),5,'d','T');
printresult(mu(:,1),5,'s','μ','см^−1');
fprintf('\nС фильтром\n');
printresult(T(:,2),5,'d','T');
printresult(mu(:,2),5,'s','μ','см^−1');
printresult(Z_K,1,'d','Z_K');
%% μ(Z)
figure('Name','μ(Z)','Units','centimeters');
hold on;

marker = {'o','s'};
color = {'b','r'};
for i = 1:numel(color)
    plot(Z,mu(:,i),marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itZ','FontName','Liberation Serif');
ylabel('\itμ\rm, см^{−1}','FontName','Liberation Serif');
legend('без фильтра','с фильтром','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 3.
unpackstruct(part3);

% λ, пм
wl = 2.*D.*sind(theta);
T = I(:,2:end)./I(:,1);
% μ, см^−1
mu = -log(T)./d.*10;
% μ_m, см^2/г
mu_m = mu./rho;
% σ_m, см^2/г
sigma_m = 8.*pi./3.*r_e.^2.*Z.*N_A./M.*1e4;
% τ_m, см^2/г
tau_m = mu_m - sigma_m;
% Аппроксимация ln(τ_m) = f(ln(λ))
[x, y, fit3] = deal(cell(size(Z)));
[a, b, d_a, d_b] = deal(zeros(size(fit3)));
for i = 1:numel(fit3)
    x{i} = log(wl(1:num_fit_points{i}));
    y{i} = log(tau_m(1:num_fit_points{i},i));
    fit3{i} = fit(x{i},y{i},{'x','1'});
    [a(i), b(i)] = matsplit(coeffvalues(fit3{i}));
    [d_a(i), d_b(i)] = matsplit(coeffvariances(fit3{i}));
end
% α
alpha = a;
d_alpha = d_a;
% k_1
k_1 = exp(b);
d_k_1 = abs(exp(b).*d_b);
%% Вывод результатов
printresult(a,d_a,'()','a');
printresult(b,d_b,'()','b');
printresult(alpha,d_alpha,'()','α');
printresult(k_1,d_k_1,'()','k_1',[],'auto');
%% I(θ)
figure('Name','I(θ)','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:numel(color)
    plot(theta,I(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itθ\rm, °','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
legend('прямой пучок','медная фольга','циркониевая фольга','Location','northwest');
xlim([min(theta) max(theta)]);
pbaspect([2 1 1]);
figresize('width',17);
%% I(λ)
figure('Name','I(λ)','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:numel(color)
    plot(wl,I(:,i),'Color',color{i});
end

set(gca,'box','on','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, пм','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
legend('прямой пучок','медная фольга','циркониевая фольга','Location','northwest');
xlim([min(wl) max(wl)]);
pbaspect([2 1 1]);
figresize('width',17);
%% T(λ)
figure('Name','T(λ)','Units','centimeters');
hold on;

color = {'r','g'};
for i = 1:numel(color)
    plot(wl,T(:,i),'Color',color{i});
end

set(gca,'box','on','FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, пм','FontName','Liberation Serif');
ylabel('\itT','FontName','Liberation Serif');
legend('медная фольга','циркониевая фольга','Location','northwest');
xlim([min(wl) max(wl)]);
pbaspect([2 1 1]);
figresize('width',17);
%% τ_m(λ)
figure('Name','τ_m(λ)','Units','centimeters');
hold on;

color = {'r','g'};
for i = 1:numel(color)
    plot(wl,tau_m(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, пм','FontName','Liberation Serif');
ylabel('\itτ_m\rm, см^2/г','FontName','Liberation Serif');
legend('медная фольга','циркониевая фольга','Location','northwest');
xlim([min(wl) max(wl)]);
pbaspect([2 1 1]);
figresize('width',17);
%% ln(τ_m) = f(ln(λ))
figure('Name','ln(τ_m) = f(ln(λ))','Units','centimeters');
hold on;

marker = {'s','^'};
color = {'r','g'};
subset = zeros(size(color));
for i = 1:numel(color)
    X = x{i}([1 end]);
    plot(X,fit3{i}(X),'Color',color{i});
    subset(i) = plot(x{i},y{i},marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('ln(\itλ\rm), ln(пм)','FontName','Liberation Serif');
ylabel('ln(\itτ_m\rm), ln(см^2/г)','FontName','Liberation Serif');
legend(subset,'медная фольга','циркониевая фольга','Location','northwest');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 4.
unpackstruct(part4);

% λ, пм
wl = 2.*D.*sind(theta);
T = I(2:end)./I(1);
% μ, см^−1
mu = -log(T)./d.*10;
% μ_m, см^2/г
mu_m = mu./rho;
% σ_m, см^2/г
sigma_m = 8.*pi./3.*r_e.^2.*Z.*N_A./M.*1e4;
% τ_m, см^2/г
tau_m = mu_m - sigma_m;
% Аппроксимация ln(τ_m) = f(ln(Z))
x = log(Z);
y = log(tau_m);
fit4 = fit(x,y,{'x','1'});
[a, b] = matsplit(coeffvalues(fit4));
[d_a, d_b] = matsplit(coeffvariances(fit4));
% β
beta = a;
d_beta = d_a;
% k_2
k_2 = exp(b);
d_k_2 = abs(exp(b).*d_b);
%% Вывод результатов
printresult(wl,2,'s','λ','пм');
printresult(T,5,'d','T');
printresult(tau_m,3,'s','τ_m','см^2/г');
printresult(a,d_a,'()','a');
printresult(b,d_b,'()','b');
printresult(beta,d_beta,'()','β');
printresult(k_2,d_k_2,'()','k_2',[],'auto');
%% τ_m(Z)
figure('Name','τ_m(Z)','Units','centimeters');
plot(Z,tau_m,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itZ','FontName','Liberation Serif');
ylabel('\itτ_m\rm, см^2/г','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% ln(τ_m) = f(ln(Z))
figure('Name','ln(τ_m) = f(ln(Z))','Units','centimeters');
hold on;

X = [2.5 4];
plot(X,fit4(X),'Color','g');
plot(x,y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('ln(\itZ\rm)','FontName','Liberation Serif');
ylabel('ln(\itτ_m\rm), ln(см^2/г)','FontName','Liberation Serif');
pbaspect([1 1 1]);
grid on;
figresize('width',9);