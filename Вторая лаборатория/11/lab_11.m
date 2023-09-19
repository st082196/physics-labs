load('data_11');

% n, °
n = (n_1 + n_2)/2;
d_n = sqrt(d_n_1^2 + d_n_2^2)/2;
% λ, Å
wl = interp1(n,wl_n,phi);
% J_фот(λ), мкА
J_photo_wl = U_photo_wl/R*1e3;
% r_λ, Вт/м^3
func_r_wl = @(wl) 8.*pi.*h.*c./(wl.*1e-10).^5.*1./(exp(h.*c./(wl.*1e-10.*k.*(T+273.15))) - 1);
r_wl = func_r_wl(wl);
% J_фот(U_пит), мкА
J_photo_supply = U_photo_supply/R*1e3;
% J_0, нА
J_0 = U_0/R*1e3;
%% Вывод результатов
printangle(n_1,d_n_1,'n_1');
printangle(n_2,d_n_2,'n_2');
printangle(n,d_n,'n');
%% λ(n)
figure('Name','λ(n)','Units','centimeters');
plot(n,wl_n,'o-','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itn\rm, °','FontName','Liberation Serif');
ylabel('\itλ\rm, Å','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% J_фот(λ)
figure('Name','J_фот(λ)','Units','centimeters');
plot(wl,J_photo_wl,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, Å','FontName','Liberation Serif');
ylabel('\itJ\rm_{фот}, мкА','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% r_λ(λ)
figure('Name','r_λ(λ)','Units','centimeters');
hold on;

X = 4000:10:7000;
plot(X,func_r_wl(X),'Color','b');
plot(wl,r_wl,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, Å','FontName','Liberation Serif');
ylabel('\itr_λ\rm, Вт/м^3','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% J_фот/r_λ = f(λ)
figure('Name','J_фот/r_λ = f(λ)','Units','centimeters');
plot(wl,J_photo_wl./r_wl,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, Å','FontName','Liberation Serif');
ylabel('\itJ\rm_{фот}/\itr_λ\rm, мкА·м^3/Вт','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% J_фот(U_пит)
figure('Name','J_фот(U_пит)','Units','centimeters');
plot(U_supply,J_photo_supply,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_{пит}, В','FontName','Liberation Serif');
ylabel('\itJ\rm_{фот}, мкА','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('height',12);
%% J_0(U_пит)
figure('Name','J_0(U_пит)','Units','centimeters');
plot(U_supply_0,J_0,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);
set(gca,'FontName','Calibri','FontSize',14,'Box','off','LooseInset',get(gca,'TightInset'));
xlabel('\itU\rm_{пит}, В','FontName','Liberation Serif');
ylabel('\itJ\rm_0, нА','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',10);