load('data_9');
%% 1.
unpackstruct(part1);

d_theta_syst = dms2degrees(d_theta_syst_DMS);
% θ_⊥, °
sample_theta_perp = dms2degrees(sample_theta_perp_DMS);
[theta_perp, d_theta_perp] = samplemean(sample_theta_perp,d_theta_syst);
% d_theta_perp = dms2degrees(d_theta_DMS);
% θ_0, °
sample_theta_0 = dms2degrees(sample_theta_0_DMS);
[theta_0, d_theta_0] = samplemean(sample_theta_0,d_theta_syst);
% ψ, °
psi = abs(theta_0 - theta_perp);
d_psi = sqrt(d_theta_0.^2 + d_theta_perp.^2);
% θ, °
theta = dms2degrees(theta_DMS);
d_theta = dms2degrees(d_theta_DMS);
% ϕ, °
phi.calc = abs(theta - theta_perp);
d_phi.calc = sqrt((d_theta).^2 + (d_theta_perp).^2);
% D_расч, ″/Å
D.calc = rad2deg((sind(psi) - sind(phi.calc))./(wl.*cosd(phi.calc))).*3600;
d_D.calc = rad2deg(1./(wl.*cosd(phi.calc)).*sqrt( ...
    (cosd(psi).*deg2rad(d_psi)).^2 + ...
    ((sind(phi.calc).*sind(psi)-1)./cosd(phi.calc).*deg2rad(d_phi.calc)).^2 + ...
    ((sind(psi)-sind(phi.calc))./wl.*d_wl).^2 ...
    )).*3600;

[k.exp, phi.exp, delta_phi, delta_wl] = deal([]);
for i = 1:max(k.calc)
    % k
    k_slice = k.calc(k.calc==i);
    k.exp = cat(1,k.exp,k_slice(1:end-1));
    % ϕ, °
    phi_slice = phi.calc(k.calc==i);
    phi.exp = cat(1,phi.exp,(phi_slice(1:end-1)+phi_slice(2:end))./2);
    % ∆ϕ, °
    delta_phi = cat(1,delta_phi,diff(phi.calc(k.calc==i)));
    % ∆λ, Å
    delta_wl = cat(1,delta_wl,diff(wl(k.calc==i)));
end
% D_эксп, ″/Å
D.exp = -3600.*delta_phi./delta_wl;
d_D.exp = sqrt(2).*abs(D.exp).*sqrt( ...
    (d_phi.calc./delta_phi).^2 + ...
    (d_wl./delta_wl).^2 ...
    );
% Аппроксимация y(x)
x = sind(psi) - sind(phi.calc);
y = k.calc.*wl*1e-4; % мкм
fit_x_y = fit(x,y,{'x'});
% d, мкм
d = coeffvalues(fit_x_y);
d_d = coeffvariances(fit_x_y);
%% Вывод результатов
printangle(theta_perp,d_theta_perp,'θ_⊥');
printangle(theta_0,d_theta_0,'θ_0');
printangle(psi,d_psi,'ψ');
printangle(theta,d_theta,'θ');
printangle(phi.calc,d_phi.calc,'ϕ');
printresult(D.calc,d_D.calc,'pm','D_расч','″/Å');
printangle(phi.exp,0,'ϕ');
printangle(delta_phi,0,'∆ϕ');
printresult(delta_wl,1,'d','∆λ','Å');
printresult(D.exp,d_D.exp,'pm','D_эксп','″/Å');
printresult(d,d_d,'()','d','мкм');
%% D(ϕ)
figure('Name','D(ϕ)','Units','centimeters');
hold on;

marker = {'o','s','^'};
color = {'b','r','g'};
subset = zeros(size(color));
for i = 1:numel(color)
    subset(i) = plot(phi.calc(k.calc==i),D.calc(k.calc==i),marker{i},'Color',color{i},'MarkerFaceColor',color{i},'MarkerSize',3,'DisplayName',"\itk\rm = "+i);
    plot(phi.exp(k.exp==i),D.exp(k.exp==i),marker{i},'Color',color{i},'MarkerSize',3);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itϕ\rm, °','FontName','Liberation Serif');
ylabel('\itD\rm, ″/Å','FontName','Liberation Serif');
legend(subset,'Location','northeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% y(x)
figure('Name','y(x)','Units','centimeters');
hold on;

X = [0.2 0.8];
plot(X,fit_x_y(X),'Color','g');
plot(x,y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itx','FontName','Liberation Serif');
ylabel('\ity\rm, мкм','FontName','Liberation Serif');
xlim(X);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 2.
unpackstruct(part2);

% R_λ
R_wl = mean(wl)./diff(wl);
d_R_wl = sqrt(wl(1).^2 + wl(2).^2)./diff(wl).^2.*d_wl;
% θ, °
theta = dms2degrees(theta_DMS);
% ϕ, °
phi = abs(theta - theta_perp);
d_phi = d_phi.calc;
% x, мм
[x, d_x] = samplemean(sample_x,d_x_syst);
% N
N = (x - x_0)./(d.*cosd(phi)).*1e3;
d_N = abs(N).*sqrt( ...
    (d_x./(x - x_0)).^2 + ...
    (d_x_0./(x - x_0)).^2 + ...
    (d_d./d).^2 + ...
    (tand(phi).*deg2rad(d_phi)).^2 ...
    );
% R_N
R_N = k.*N;
d_R_N = abs(k).*d_N;
%% Вывод результатов
printresult(R_wl,d_R_wl,'pm','R_λ');
printangle(theta,d_theta,'θ');
printangle(phi,d_phi,'ϕ');
printresult(x,d_x,'pm','x','мм');
printresult(N,d_N,'pm','N');
printresult(R_N,d_R_N,'pm','R_N');