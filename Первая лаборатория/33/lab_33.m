load('data_33');
%% 1.
unpackstruct(part1);

func_B = @(r) (mu_0.*I)./(2.*pi.*r).*1e6;
B.theor = func_B(r); % мТл
%% Вывод результатов
printresult(B.theor,3,'d','B_теор','мТл');
%% B(r)
figure('Name','B(r)','Units','centimeters');
hold on;

X = (0:0.1:25).';
plot(X,func_B(X),'Color','g');
plot(r,B.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itr\rm, мм','FontName','Liberation Serif');
ylabel('\itB\rm, мТл','FontName','Liberation Serif');
ylim([0 0.8]);
pbaspect([4 3 1]);
grid on;
figresize('width',10.5);
%% 2.
unpackstruct(part2);

func_B = @(z) (mu_0.*I.*R.^2)./(2.*(z.^2 + R.^2).^(3./2)).*1e5;
B.theor = func_B(z); % мТл
%% Вывод результатов
printresult(B.theor,3,'d','B_теор','мТл');
%% B(z)
figure('Name','B(z)','Units','centimeters');
hold on;

X = (-5:0.02:5).';
plot(X,func_B(X),'Color','g');
plot(z,B.exp,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itz\rm, см','FontName','Liberation Serif');
ylabel('\itB\rm, мТл','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',10.5);
%% 3.
unpackstruct(part3);

func_B = @(z) mu_0.*I.*N.*( ...
    1./sqrt((2.*z+L).^2 + 4.*(R.*1e-1).^2).*(z./L + 1./2) - ...
    1./sqrt((2.*z-L).^2 + 4.*(R.*1e-1).^2).*(z./L - 1./2) ...
    ).*1e5;
B.theor = func_B(z); % мТл
%% Вывод результатов
for j = 1:numel(L)
    fprintf('\nL = %.0f см',L(j));
    printresult(B.theor(:,j),2,'d','B_теор','мТл');
end
%% B(z)
figure('Name','B(z)','Units','centimeters');
hold on;

X = (-15:0.1:15).';
Y = func_B(X);
marker = {'o','s','^','v'};
color = {'b','r','g','m'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(X,Y(:,j),'Color',color{j});
    subset(j) = plot(z,B.exp(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itL\\rm = %.0f см',L(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itz\rm, см','FontName','Liberation Serif');
ylabel('\itB\rm, мТл','FontName','Liberation Serif');
legend(subset,'Location','northeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% 4.
unpackstruct(part4);

func_B = @(z) mu_0.*I.*N.*R.^2./2.*( ...
    1./((z.*1e1-R./2).^2 + R.^2).^(3./2) + ...
    1./((z.*1e1+R./2).^2 + R.^2).^(3./2) ...
    ).*1e6;
B.theor = func_B(z); % мТл
%% Вывод результатов
for j = 1:numel(I)
    fprintf('\nI = %.1f А',I(j));
    printresult(B.theor(:,j),2,'d','B_теор','мТл');
end
%% B(z)
figure('Name','B(z)','Units','centimeters');
hold on;

X = (-15:0.1:15).';
Y = func_B(X);
marker = {'o','s','^'};
color = {'b','r','g'};
subset = zeros(size(color));
for j = 1:numel(color)
    plot(X,Y(:,j),'Color',color{j});
    subset(j) = plot(z,B.exp(:,j),marker{j},'Color',color{j},'MarkerFaceColor',color{j},'MarkerSize',3,'DisplayName',sprintf('\\itI\\rm = %.1f А',I(j)));
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itz\rm, см','FontName','Liberation Serif');
ylabel('\itB\rm, мТл','FontName','Liberation Serif');
legend(subset,'Location','northeast','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);