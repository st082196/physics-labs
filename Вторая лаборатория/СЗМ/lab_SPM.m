load('data_SPM')

% z, нм
z = reshape(data,[],1)*0.54;
% z_pit, нм
sample_z_pit = z(z < 30);
[z_pit, d_z_pit] = samplemean(sample_z_pit,d_z_syst);
% z_surf, нм
sample_z_surf = z(z > 120);
[z_surf, d_z_surf] = samplemean(sample_z_surf,d_z_syst);
% h, нм
h = z_surf - z_pit;
d_h = sqrt(d_z_pit.^2 + d_z_surf.^2);
% l, нм
[l, d_l] = samplemean(sample_l,d_l_syst);
% w, нм
[w, d_w] = samplemean(sample_w,d_w_syst);
% ρ, бит/мм^2
rho = 1.5./(l.*x).*1e12;
d_rho = abs(rho).*sqrt((d_l./l).^2 + (d_x./x).^2);
% V, ГБайт
V = rho.*S./8.*1e-9;
d_V = abs(S.*d_rho./8.*1e-9);
%% Вывод результатов
printresult(z_pit,d_z_pit,'pm','z_pit','нм');
printresult(z_surf,d_z_surf,'pm','z_surf','нм');
printresult(h,d_h,'pm','h','нм');
printresult(l,d_l,'pm','l','нм');
printresult(w,d_w,'pm','w','нм');
printresult(rho,d_rho,'pm','ρ','бит/мм^2','auto');
printresult(V,d_V,'pm','V','Гбайт');
%% z(№ пикселя)
figure('Name','z(№ пикселя)','Units','centimeters');

plot(z,'.','Color','#B39719');
yline(30,'--','\itz\rm = 30','FontName','Calibri');
yline(120,'--','\itz\rm = 120','FontName','Calibri');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('№ пикселя');
ylabel('\itz','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);