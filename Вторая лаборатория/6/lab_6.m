load('data_6');

% E, В/м
E = U./d.*1e2;
% δ, °
delta = 2.*(phi - phi(1,:));
% Аппроксимация  δ(E^2)
fit1 = fit(repmat(E.^2,size(delta,2),1),reshape(delta,[],1),{'x'});
a = coeffvalues(fit1); % м^2⋅°/В^2
d_a = coeffvariances(fit1);
% B, м/В^2
B = a./(360.*l).*1e2;
d_B = d_a./(360.*l).*1e2;
%% Вывод результатов
for i = 1:size(delta,2)
    printresult(delta(:,i),2,'d',"δ (серия №"+i+")",'°');
end
printresult(a,d_a,'pm','a','м^2⋅°/В^2','auto');
printresult(B,d_B,'pm','B','м/В^2','auto');
%% δ(E^2)
figure('Name','δ(E^2)','Units','centimeters');
hold on;

X = [0 4e11];
plot(X,fit1(X),'Color','g','DisplayName','аппроксимирующая прямая');
marker = {'o','s','^'};
for i = 1:numel(marker)
    plot(E.^2,delta(:,i),marker{i},'Color','b','MarkerFaceColor','b','MarkerSize',3,'DisplayName',"серия №"+i);
end

set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
xlabel('\itE\rm^2, В/м','FontName','Liberation Serif');
ylabel('\itδ\rm, °','FontName','Liberation Serif');
legend('Location','southeast');
pbaspect([4 3 1]);
grid on;
figresize('width',17);