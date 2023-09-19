load('data_29');
circuit = 'RCserial';

switch circuit
    case 'RCserial'
        unpackstruct(RCserial);
        element = 'C';
        func_r_f = @(f) 1./(2.*pi.*f.*C.*1e-6);
    case 'RLserial'
        unpackstruct(RLserial);
        element = 'L';
        func_r_f = @(f) 2.*pi.*f.*L.*1e-3;
end

% φ(f,∆t), °
phi_f_delta_t = f.*delta_t.*1e-3.*360;
d_phi_f_delta_t = sqrt((delta_t.*d_f).^2 + (f.*d_delta_t).^2).*1e-3.*360;
% φ(U_0,I_0), °
phi_U_0_I_0 = sign(delta_t(1)).*acosd(I_0.*r./U_0.*1e-3);
d_phi_U_0_I_0 = rad2deg(((U_0./(I_0.*1e-3.*r)).^2 - 1).^-0.5.*sqrt((d_U_0./U_0).^2 + (d_I_0./I_0).^2));
% φ(f,C/L), °
func_phi_f = @(f) sign(delta_t(1)).*atand(func_r_f(f)./r);
phi_f = func_phi_f(f);

% |Z|(φ), Ом
Z_phi = r./cosd(phi_f_delta_t);
d_Z_phi = abs(r.*sind(phi_f_delta_t)./cosd(phi_f_delta_t).^2.*deg2rad(d_phi_f_delta_t));
% |Z|(U_0,I_0), Ом
Z_U_0_I_0 = U_0./I_0.*1e3;
d_Z_U_0_I_0 = abs(Z_U_0_I_0).*sqrt((d_U_0./U_0).^2 + (d_I_0./I_0).^2);
% |Z|(f,C/L), Ом
func_Z_f = @(f) sqrt(r.^2 + func_r_f(f).^2);
Z_f = func_Z_f(f);

% r_C/L(φ), Ом
r_phi = sign(delta_t(1)).*r.*tand(phi_f_delta_t);
d_r_phi = abs(r./cosd(phi_f_delta_t).^2.*deg2rad(d_phi_f_delta_t));
% r_C/L(U_0,I_0), Ом
r_U_0_I_0 = sqrt((U_0./I_0.*1e3).^2 - r.^2);
d_r_U_0_I_0 = (U_0./I_0.*1e3).^2.*((U_0./I_0.*1e3).^2 - r.^2).^-0.5.*sqrt((d_U_0./U_0).^2 + (d_I_0./I_0).^2);
% r_C/L(f,C/L), Ом
r_f = func_r_f(f);

% Аппроксимация y(x)
switch circuit
    case 'RCserial'
        x = tand(-phi_f_delta_t);
        y = 1./(2.*pi.*f.*r).*1e6; % мкФ
        fit_x_y.f = fit(x,y,{'x'});
        fit_x_y.C = coeffvalues(fit_x_y.f);
        fit_x_y.d_C = coeffvariances(fit_x_y.f);
    case 'RLserial'
        x = 2.*pi.*f./r.*1e-3; % мГн^−1
        y = tand(phi_f_delta_t);
        fit_x_y.f = fit(x,y,{'x'});
        fit_x_y.L = coeffvalues(fit_x_y.f);
        fit_x_y.d_L = coeffvariances(fit_x_y.f);
end
%% Вывод результатов
printresult(phi_f_delta_t,d_phi_f_delta_t,'pm','φ(f,∆t)','°');
printresult(phi_U_0_I_0,d_phi_U_0_I_0,'pm','φ(U_0,I_0)','°');
printresult(phi_f,3,'s',['φ(f,',element,')'],'°');

printresult(Z_phi,d_Z_phi,'pm','|Z|(φ)','Ом');
printresult(Z_U_0_I_0,d_Z_U_0_I_0,'pm','|Z|(U_0,I_0)','Ом');
printresult(Z_f,3,'s',['|Z|(f,',element,')'],'Ом');

printresult(r_phi,d_r_phi,'pm',['r_',element,'(φ)'],'Ом');
printresult(r_U_0_I_0,d_r_U_0_I_0,'pm',['r_',element,'(U_0,I_0)'],'Ом');
printresult(r_f,3,'s',['r_',element,'(f,',element,')'],'Ом');

switch circuit
    case 'RCserial'
        printresult(fit_x_y.C,fit_x_y.d_C,'pm','C','мкФ');
    case 'RLserial'
        printresult(fit_x_y.L,fit_x_y.d_L,'pm','L','мГн');
end
%% φ(f)
figure('Name','φ(f)','Units','centimeters');
hold on;

X = 0:10:1000;
plot(X,func_phi_f(X),'Color','g','DisplayName',['\itφ\rm(\itf\rm,\it',element,'\rm)']);
errorbar(f,phi_U_0_I_0,d_phi_U_0_I_0,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'CapSize',3,'DisplayName','\itφ\rm(\itU\rm_0,\itI\rm_0)');
errorbar(f,phi_f_delta_t,d_phi_f_delta_t,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'CapSize',3,'DisplayName','\itφ\rm(\itf\rm,\it∆t\rm)');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\itφ\rm, °','FontName','Liberation Serif');
legend('Location','southeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% |Z|(f)
figure('Name','|Z|(f)','Units','centimeters');
hold on;

X = 0:10:1000;
plot(X,func_Z_f(X),'Color','g','DisplayName',['|\itZ\rm|(\itf\rm,\it',element,'\rm)']);
errorbar(f,Z_U_0_I_0,d_Z_U_0_I_0,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'CapSize',3,'DisplayName','|\itZ\rm|(\itU\rm_0,\itI\rm_0)');
errorbar(f,Z_phi,d_Z_phi,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'CapSize',3,'DisplayName','|\itZ\rm|(\itφ\rm)');

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel('\it|Z|\rm, Ом','FontName','Liberation Serif');
legend('Location','northeast','FontName','Liberation Serif');
ylim([100 400]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% r_C/L(f)
figure('Name',['r_',element,'(f)'],'Units','centimeters');
hold on;

X = 0:10:1000;
plot(X,func_r_f(X),'Color','g','DisplayName',['\itr_',element,'\rm(\itf\rm,\it',element,'\rm)']);
errorbar(f,r_U_0_I_0,d_r_U_0_I_0,'s','Color','r','MarkerFaceColor','r','MarkerSize',3,'CapSize',3,'DisplayName',['\itr_',element,'\rm(\itU\rm_0,\itI\rm_0)']);
errorbar(f,r_phi,d_r_phi,'o','Color','b','MarkerFaceColor','b','MarkerSize',3,'CapSize',3,'DisplayName',['\itr_',element,'\rm(\itφ\rm)']);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itf\rm, Гц','FontName','Liberation Serif');
ylabel(['\itr_',element,'\rm, Ом'],'FontName','Liberation Serif');
legend('Location','northeast','FontName','Liberation Serif');
ylim([0 400]);
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% y(x)
figure('Name','y(x)','Units','centimeters');
hold on;

X = [0 3.5];
plot(X,fit_x_y.f(X),'Color','g');
plot(x,y,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
switch circuit
    case 'RCserial'
        xlabel('\itx','FontName','Liberation Serif');
        ylabel('\ity\rm, мкФ','FontName','Liberation Serif');
    case 'RLserial'
        xlabel('\itx\rm, мГн^{−1}','FontName','Liberation Serif');
        ylabel('\ity','FontName','Liberation Serif');
end
pbaspect([4 3 1]);
grid on;
figresize('width',17);