load('data_18'); %#ok<*SAGROW>

% Линзы по отдельности
for j = 1:numel(lens)
    % u, см
    lens(j).u = lens(j).x - lens(j).q;
    % b, см
    lens(j).b = lens(j).r - lens(j).q;
    % γ
    lens(j).gamma = lens(j).y_img./lens(j).y;
    % Аппроксимация γ(b)
    lens(j).fit1.object = fit(lens(j).b,lens(j).gamma,{'x','1'});
    [lens(j).fit1.a, lens(j).fit1.b] = matsplit(coeffvalues(lens(j).fit1.object));
    [lens(j).fit1.d_a, lens(j).fit1.d_b] = matsplit(coeffvariances(lens(j).fit1.object));
    % p_з, дптр
    lens(j).p_rear = -lens(j).fit1.a.*1e2;
    lens(j).d_p_rear = lens(j).fit1.d_a.*1e2;
    % Аппроксимация γ^−1(u)
    lens(j).fit2.object = fit(lens(j).u,lens(j).gamma.^-1,{'x','1'});
    [lens(j).fit2.a, lens(j).fit2.b] = matsplit(coeffvalues(lens(j).fit2.object));
    [lens(j).fit2.d_a, lens(j).fit2.d_b] = matsplit(coeffvariances(lens(j).fit2.object));
    % p_п, дптр
    lens(j).p_front = lens(j).fit2.a.*1e2;
    lens(j).d_p_front = lens(j).fit2.d_a.*1e2;
    % p, дптр
    lens(j).p = (lens(j).p_front + lens(j).p_rear)./2;
    lens(j).d_p = 1./2.*sqrt(lens(j).d_p_front.^2 + lens(j).d_p_rear.^2);
    % f, см
    lens(j).f = 1./lens(j).p.*1e2;
    lens(j).d_f = lens(j).d_p./lens(j).p.^2.*1e2;
end

% L, мм
L = d - l(1) - l(2);
d_L = sqrt(d_d.^2 + d_l(1).^2 + d_l(2).^2);
% A (теор.)
A.theor = 1 - lens(1).p.*L.*1e-3;
d_A.theor = sqrt((L.*lens(1).d_p).^2 + (lens(1).p.*d_L).^2).*1e-3;
% B (теор.), мм
B.theor = L;
d_B.theor = d_L;
% C (теор.), дптр
C.theor = -(lens(1).p + lens(2).p - lens(1).p.*lens(2).p.*L.*1e-3);
d_C.theor = sqrt( ...
    ((lens(2).p.*L.*1e-3 - 1).*lens(1).d_p).^2 + ...
    ((lens(1).p.*L.*1e-3 - 1).*lens(2).d_p).^2 + ...
    (lens(1).p.*lens(2).p.*d_L.*1e-3).^2 ...
    );
% D (теор.)
D.theor = 1 - lens(2).p.*L.*1e-3;
d_D.theor = sqrt((L.*lens(2).d_p).^2 + (lens(2).p.*d_L).^2).*1e-3;
% p (теор.), дптр
p.theor = C.theor;
d_p.theor = d_C.theor;
% f (теор.), см
f.theor = 1./p.theor.*1e2;
d_f.theor = d_p.theor./p.theor.^2.*1e2;

% Оптическая система
for j = 1:numel(system)
    % u, см
    system(j).u = system(j).x - system(j).q;
    % b, см
    system(j).b = system(j).r - system(j).q;
    % γ
    system(j).gamma = system(j).y_img./system(j).y;
    % Аппроксимация γ(b)
    system(j).fit1.object = fit(system(j).b,system(j).gamma,{'x','1'});
    [system(j).fit1.a, system(j).fit1.b] = matsplit(coeffvalues(system(j).fit1.object));
    [system(j).fit1.d_a, system(j).fit1.d_b] = matsplit(coeffvariances(system(j).fit1.object));
    % C_1, дптр
    system(j).C_1 = system(j).fit1.a.*1e2;
    system(j).d_C_1 = system(j).fit1.d_a.*1e2;
    % A
    system(j).A = system(j).fit1.b;
    system(j).d_A = system(j).fit1.d_b;
    % Аппроксимация γ^−1(u)
    system(j).fit2.object = fit(system(j).u,system(j).gamma.^-1,{'x','1'});
    [system(j).fit2.a, system(j).fit2.b] = matsplit(coeffvalues(system(j).fit2.object));
    [system(j).fit2.d_a, system(j).fit2.d_b] = matsplit(coeffvariances(system(j).fit2.object));
    % C_2, дптр
    system(j).C_2 = -system(j).fit2.a.*1e2;
    system(j).d_C_2 = system(j).fit2.d_a.*1e2;
    % E
    system(j).E = system(j).fit2.b;
    system(j).d_E = system(j).fit2.d_b;
end

% A (эксп.)
A.exp = system(1).A;
d_A.exp = system(1).d_A;
% D (эксп.)
D.exp = system(2).A;
d_D.exp = system(2).d_A;
% C (эксп.), дптр
C.exp = (system(1).C_1 + system(1).C_2 + system(2).C_1 + system(2).C_2)./4;
d_C.exp = 1./4.*sqrt(system(1).d_C_1.^2 + system(1).d_C_2.^2 + system(2).d_C_1.^2 + system(2).d_C_2.^2);
% B (эксп.), мм
B.exp = (system(1).A + system(2).A - system(1).E - system(2).E)./(2.*C.exp).*1e3;
d_B.exp = abs(1./(2.*C.exp)).*sqrt( ...
    system(1).d_A.^2 + system(2).d_A.^2 + system(1).d_E.^2 + system(2).d_E.^2 + ...
    ((system(1).A + system(2).A - system(1).E - system(2).E)./C.exp.*d_C.exp).^2 ...
    ).*1e3;
% p (эксп.), дптр
p.exp = C.exp;
d_p.exp = d_C.exp;
% f (эксп.), см
f.exp = 1./p.exp.*1e2;
d_f.exp = d_p.exp./p.exp.^2.*1e2;
%% Вывод результатов
for j = 1:numel(lens)
    fprintf('\n\nЛинза №%u\n',j);
    printresult(lens(j).u,2,'d','u','см');
    printresult(lens(j).b,2,'d','b','см');
    printresult(lens(j).gamma,3,'d','γ');
    fprintf('\nАппроксимация γ(b)\n');
    printresult(lens(j).fit1.a,lens(j).fit1.d_a,'pm','a','см^−1');
    printresult(lens(j).fit1.b,lens(j).fit1.d_b,'pm','b');
    printresult(lens(j).p_rear,lens(j).d_p_rear,'pm','p_з','дптр');
    fprintf('\nАппроксимация γ^−1(u)\n');
    printresult(lens(j).fit2.a,lens(j).fit2.d_a,'pm','a','см^−1');
    printresult(lens(j).fit2.b,lens(j).fit2.d_b,'pm','b');
    printresult(lens(j).p_front,lens(j).d_p_front,'pm','p_п','дптр');
    printresult(lens(j).p,lens(j).d_p,'pm','p','дптр');
    printresult(lens(j).f,lens(j).d_f,'pm','f','см');
end

fprintf('\n\nТеоретический расчёт оптической системы\n');
printresult(L,d_L,'pm','L','мм');
printresult(A.theor,d_A.theor,'()','A');
printresult(B.theor,d_B.theor,'()','B','мм');
printresult(C.theor,d_C.theor,'()','C','дптр');
printresult(D.theor,d_D.theor,'()','D');
printresult(p.theor,d_p.theor,'pm','p','дптр');
printresult(f.theor,d_f.theor,'pm','f','см');

for j = 1:numel(system)
    fprintf('\n\nОптическая система №%u\n',j);
    printresult(system(j).u,1,'d','u','см');
    printresult(system(j).b,1,'d','b','см');
    printresult(system(j).gamma,3,'d','γ');
    fprintf('\nАппроксимация γ(b)\n');
    printresult(system(j).fit1.a,system(j).fit1.d_a,'pm','a','см^−1');
    printresult(system(j).fit1.b,system(j).fit1.d_b,'pm','b');
    printresult(system(j).C_1,system(j).d_C_1,'pm','C_1','дптр');
    printresult(system(j).A,system(j).d_A,'pm','A');
    fprintf('\nАппроксимация γ^−1(u)\n');
    printresult(system(j).fit2.a,system(j).fit2.d_a,'pm','a','см^−1');
    printresult(system(j).fit2.b,system(j).fit2.d_b,'pm','b');
    printresult(system(j).C_2,system(j).d_C_2,'pm','C_2','дптр');
    printresult(system(j).E,system(j).d_E,'pm','E');
end

fprintf('\n\nЭкспериментальный расчёт оптической системы\n');
printresult(A.exp,d_A.exp,'()','A');
printresult(B.exp,d_B.exp,'()','B','мм');
printresult(C.exp,d_C.exp,'()','C','дптр');
printresult(D.exp,d_D.exp,'()','D');
printresult(p.exp,d_p.exp,'pm','p','дптр');
printresult(f.exp,d_f.exp,'pm','f','см');
%% γ(b)
unpackstruct(lens(1));

figure('Name','γ(b)','Units','centimeters');
hold on;

x_fit = [16 32];
plot(x_fit,fit1.object(x_fit),'Color','g');
plot(b,gamma,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itb\rm, см','FontName','Liberation Serif');
ylabel('\itγ','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);
%% γ^−1(u)
unpackstruct(lens(1));

figure('Name','γ^−1(u)','Units','centimeters');
hold on;

x_fit = [-42 -30];
plot(x_fit,fit2.object(x_fit),'Color','g');
plot(u,gamma.^-1,'o','Color','b','MarkerFaceColor','b','MarkerSize',3);

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itu\rm, см','FontName','Liberation Serif');
ylabel('\itγ\rm^{−1}','FontName','Liberation Serif');
pbaspect([4 3 1]);
grid on;
figresize('width',17);