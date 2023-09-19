load('data_46');
%% 1.
unpackstruct(part1);

I = samplemean(sample_I); % мА
I_0 = max(I); % мА
phi_0 = phi(I == I_0); % °
f.exp = I./I_0;
f.theor = cosd(phi - phi_0).^2;
%% Вывод результатов
printresult(I,2,'d','I','мА');
printresult(f.exp,2,'d','f_эксп');
printresult(f.theor,2,'d','f_теор');
%% f(φ − φ_0)
figure('Name','f(φ − φ_0)','Units','centimeters');
hold on;

x = linspace(-90,90,500);
y = cosd(x).^2;
plot(x,y,'Color','g');
plot(phi-phi_0,f.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itφ\rm − \itφ\rm_0, °','FontName','Liberation Serif');
ylabel('\itf','FontName','Liberation Serif');
legend('теоретическая кривая','экспериментальные точки','Location','south');
xlim(x([1 end]));
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% φ_0, °
[phi_0, d_phi_0] = samplemean(sample_phi_0,d_phi_syst);
%% Вывод результатов
printresult(phi_0,d_phi_0,'pm','φ_0','°');
%% 3.
unpackstruct(part3);

% φ, °
[phi, d_phi] = samplemean(sample_phi,d_phi_syst);
% φ′, °
phi_rot = phi - phi_0;
d_phi_rot = sqrt(d_phi.^2 + d_phi_0.^2);
% K
K.exp = phi_rot(:,2)./phi_rot(:,1);
d_K.exp = abs(K.exp).*sqrt((d_phi_rot(:,1)./phi_rot(:,1)).^2 + (d_phi_rot(:,2)./phi_rot(:,2)).^2);
K.theor = -d(2)./d(1);
%% Вывод результатов
printresult(phi(:,1),d_phi(:,1),'pm','φ_л','°');
printresult(phi(:,2),d_phi(:,2),'pm','φ_пр','°');
printresult(phi_rot(:,1),d_phi_rot(:,1),'pm','φ′_л','°');
printresult(phi_rot(:,2),d_phi_rot(:,2),'pm','φ′_пр','°');
printresult(K.exp,d_K.exp,'pm','K_эксп');
printresult(K.theor,3,'s','K_теор');
%% 4.
unpackstruct(part4);

% φ, °
[phi, d_phi] = samplemean(sample_phi,d_phi_syst);
% φ′, °
phi_rot = phi - phi_0;
d_phi_rot = sqrt(d_phi.^2 + d_phi_0.^2);
% K
K.exp = phi_rot(:,1,:)./phi_rot(:,2,:);
d_K.exp = abs(K.exp).*sqrt((d_phi_rot(:,1,:)./phi_rot(:,1,:)).^2 + (d_phi_rot(:,2,:)./phi_rot(:,2,:)).^2);
K.theor = -d(1)./d(2);
d_K.theor = abs(K.theor).*sqrt((d_d./d(1)).^2 + (d_d./d(2)).^2);
% Аппроксимация φ′(d⋅C)
x = reshape(d.*C,[],numel(d)+numel(C)).';
y = num2cell(reshape(phi_rot,[],numel(d)+numel(C)).',1).';
fit_x_y = cell(size(y));
[alpha, d_alpha] = deal(zeros(size(fit_x_y)));
for i = 1:numel(fit_x_y)
    fit_x_y{i} = fit(x,y{i},{'x'});
    % [α], °⋅см^2⋅г^−1
    alpha(i) = coeffvalues(fit_x_y{i});
    d_alpha(i) = coeffvariances(fit_x_y{i});
end
%% Вывод результатов
for k = 1:numel(C)
    printresult(C(k),2,'d','C','г/мл');
    printresult(phi(:,1,k),d_phi(:,1,k),'pm','φ_1','°');
    printresult(phi(:,2,k),d_phi(:,2,k),'pm','φ_2','°');
    printresult(phi_rot(:,1,k),d_phi_rot(:,1,k),'pm','φ′_1','°');
    printresult(phi_rot(:,2,k),d_phi_rot(:,2,k),'pm','φ′_2','°');
    printresult(K.exp(:,:,k),d_K.exp(:,:,k),'pm','K_эксп');
end
printresult(K.theor,d_K.theor,'pm','K_теор');
printresult(alpha,d_alpha,'pm','[α]','°⋅см^2⋅г^−1');
%% φ′(d⋅C)
figure('Name','φ′(d⋅C)','Units','centimeters');
hold on;

x_fit = [0 5];
marker = ["o";"s";"^"];
color = ["g";"y";"r"];
subset = zeros(size(color));
for i = 1:numel(color)
    plot(x_fit,fit_x_y{i}(x_fit),'Color',color(i));
    subset(i) = plot(x,y{i},marker(i),'Color',color(i),'MarkerFaceColor',color(i),'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itd · C\rm, г/см^2','FontName','Liberation Serif');
ylabel('\itφ′\rm, °','FontName','Liberation Serif');
legend(subset,'зелёный фильтр','жёлтый фильтр','красный фильтр','Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);