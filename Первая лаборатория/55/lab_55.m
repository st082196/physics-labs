%% Интерполяция таблицы
% eds_t_table=zeros(61,1);
T = (-200:10:400).';
[eds_T_table_fit] = polyfit(T,eds_t_table,3);
p = 775;
T_N2=-195.8+0.0109*(p-760);
T_H2O=100+0.0367*(p-760)-0.000023*(p-760)^2;
%% Расчёт ε
rho = 2.8;
l_x = [791 713 708].';
R_2 = [505.0 235.6 336.8].';
L = 1;
eds_n = 1.01795;
eds_t=rho*l_x./(R_2.*L)*eds_n;
eds_t(3)=-eds_t(3);
[eds_T_fit] = polyfit([T_H2O T_N2 0],[eds_t(1) eds_t(3) 0],2);
%% Расчёт ε №2
rho = 2.8;
l_x = [791 713 708].';
R_2 = [505.0 235.6 336.8].';
L = 1;
eds_n = 1.01795;
eds_t=rho*l_x./(R_2.*L)*eds_n;
eds_t(3)=-eds_t(3);
[eds_T_fit] = polyfit([T_H2O 231.85 T_N2],[eds_t(1) eds_t(2) eds_t(3)],2);
%% ε(T) table
figure('Name','ε(T)');
title('ε(T)');
hold on;

plot(T,eds_t_table,'o','Color','b','MarkerFaceColor','b','MarkerSize',5);

x=linspace(-200,400,1000);
plot(x,polyval(eds_T_table_fit,x),'Color','g','LineWidth',1);

xlabel('T, °C');
ylabel('ε, мВ');
pbaspect([4 3 1]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
grid on;
%% Δε(T)
figure('Name','Δε(T)');
title('Δε(T)');
hold on;
x=[T_H2O 231.85 T_N2];
d_eds=eds_t(1:3)-polyval(eds_T_table_fit,x).';
plot(x,d_eds,'o','Color','b','MarkerFaceColor','b','MarkerSize',5);
plot(0,0,'o','Color','k','MarkerSize',5);

x=linspace(-200,400,1000);
d_eds_T_fit = polyfit([T_H2O T_N2 0],[d_eds(1) d_eds(3) 0],2);
plot(x,polyval(d_eds_T_fit,x),'Color','g','LineWidth',1);

xlabel('T, °C');
ylabel('Δε, мВ');
pbaspect([4 3 1]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
grid on;
%% ε(T)
figure('Name','ε(T)');
title('ε(T)');
hold on;
plot([T_H2O 231.85 T_N2],eds_t,'o','Color','b','MarkerFaceColor','b','MarkerSize',5);

x=linspace(-200,400,1000);
plot(x,polyval(eds_T_fit,x).','Color','g','LineWidth',1);

xlabel('T, °C');
ylabel('ε, мВ');
pbaspect([4 3 1]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
grid on;
%% ε(T) & ε(T)_table (ε(T) interpolated)
figure('Name','ε(T)');
title('ε(T)');
hold on;
plot([T_H2O 231.85 T_N2],eds_t,'o','Color','b','MarkerFaceColor','b','MarkerSize',5);
plot(T,eds_t_table,'o','Color','r','MarkerFaceColor','r','MarkerSize',5);

x=linspace(-200,400,1000);
plot(x,polyval(eds_T_fit,x),'Color','b','LineWidth',1);
plot(x,polyval(eds_T_table_fit,x),'Color','r','LineWidth',1);

xlabel('T, °C');
ylabel('ε, мВ');
pbaspect([4 3 1]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
grid on;
%% ε(T) & ε(T)_table (Δε(T) interpolated)
figure('Name','ε(T)');
title('ε(T)');
hold on;
plot([T_H2O 231.85 T_N2],eds_t,'o','Color','b','MarkerFaceColor','b','MarkerSize',5);
% plot(T,eds_t_table,'s','Color','r','MarkerFaceColor','r','MarkerSize',3);

x=linspace(-200,400,1000);
plot(x,polyval(eds_T_table_fit,x)+polyval(d_eds_T_fit,x),'Color','b','LineWidth',1);
plot(x,polyval(eds_T_table_fit,x),'Color','r','LineWidth',1);

xlabel('T, °C');
ylabel('ε, мВ');
legend('Экспериментальные точки','Экспериментальная интерполирующая кривая','Табличная интерполирующая кривая','Location','southeast');
pbaspect([4 3 1]);
set(findall(gcf,'-property','FontSize'),'FontSize',20)
grid on;