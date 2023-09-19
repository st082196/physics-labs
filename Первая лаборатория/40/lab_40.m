load('data_40');

% α_1, °
sample_alpha_1 = mean(sample_alpha_1,3);
[alpha_1, d_alpha_1] = samplemean(sample_alpha_1,d_alpha_syst);
% α_2, °
sample_alpha_2 = mean(sample_alpha_2,3);
[alpha_2, d_alpha_2] = samplemean(sample_alpha_2,d_alpha_syst);
% θ, °
theta = abs(alpha_2 - alpha_1)./2;
d_theta = sqrt(d_alpha_1.^2 + d_alpha_2.^2)./2;
% φ_0, °
phi_0 = mean(phi_0,3);
% φ, °
sample_phi = mean(sample_phi,3);
[phi, d_phi] = samplemean(sample_phi,d_phi_syst);
% φ_min, °
phi_min = abs(phi - phi_0);
d_phi_min = sqrt(d_phi.^2 + d_phi_0.^2);
% n
n = sind((theta + phi_min)./2)./sind(theta./2);
d_n = abs(cosd((theta + phi_min)./2)./(2.*sind(theta./2))) .* ...
    deg2rad(sqrt(((1 - tand((theta + phi_min)./2)./tand(theta./2)).*d_theta).^2 + d_phi_min.^2));
%% Вывод результатов
printresult(alpha_1,d_alpha_1,'pm','α_1','°');
printresult(alpha_2,d_alpha_2,'pm','α_2','°');
printresult(theta,d_theta,'pm','θ','°');
printresult(phi_0,d_phi_0,'pm','φ_0','°');
printresult(phi,d_phi,'pm','φ','°');
printresult(phi_min,d_phi_min,'pm','φ_min','°');
printresult(n,d_n,'pm','n');
%% n(λ)
figure('Name','n(λ)','Units','centimeters');
plot(wl,n,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itn','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);