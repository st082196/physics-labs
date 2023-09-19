%% Обрезка изображения
img = imread('img\image.png');
img = img(910:1090,2:end-1,:);
img(:,:,1:2) = 0;
img(repmat(img(:,:,3)<30,1,1,3)) = 255;
imwrite(img,'img\image_cut_matlab.png');
%% I(x)
img = imread('img\line_1.png');
I = double(img(:,:,3));
I(I==255) = NaN;
I = mean(I,'omitnan');

figure('Name','I(x)','Units','centimeters');
plot(I,'Color','b');
set(gca,'FontName','Calibri','FontSize',14,'LooseInset',get(gca,'TightInset'));
title('Линия №1');
xlabel('\itx\rm, пикс','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
xlim([0 numel(I)]);
pbaspect([3 1 1]);
figresize('width',17);
%% Расчёты
load('data_86');

% d, пикс
d = squeeze(x(:,2,:) - x(:,1,:));
% d′, пикс
d_prime = squeeze(x_prime(:,2,:) - x_prime(:,1,:));
% dν, м^−1
[delta_nu, d_delta_nu] = deal(zeros(size(d)-[1 0]));
for k = 1:size(delta_nu,1)
    delta_nu(k,:) = 1./(2.*t.*1e-3) .* (d(k,:).^2-d_prime(k,:).^2) ./ (d(k,:).^2-d(k+1,:).^2);
    d_delta_nu(k,:) = abs(d_d./(t.*1e-3)).*sqrt( ...
        ((d(k,:).*(d_prime(k,:).^2-d(k+1,:)))./(d(k,:).^2-d(k+1,:).^2).^2).^2 + ...
        (d_prime(k,:)./(d(k,:).^2-d(k+1,:).^2)).^2 + ...
        ((d(k+1,:).*(d(k,:).^2-d_prime(k,:).^2))./(d(k,:).^2-d(k+1,:).^2).^2).^2 ...
        );
end
% mean_dν, м^−1
mean_delta_nu = mean(delta_nu);
d_mean_delta_nu = sqrt(sum(d_delta_nu.^2))./size(delta_nu,1);
% A, Дж
A = mean_delta_nu.*h.*c./(F + 1);
d_A = abs(A.*d_mean_delta_nu./mean_delta_nu);
% g(I)
g = 1836.*A.*(3.*n.^3)./(8.*h.*c.*R.*alpha_sqr.*Z.*Z_a.^2).*1e-2;
d_g = abs(g.*d_A./A);
% μ_I/μ_яд
mu_I_div_mu_nucl = I.*g;
d_mu_I_div_mu_nucl = I.*d_g;
%% Вывод результатов
for i = 1:numel(mu_I_div_mu_nucl)
    fprintf('\n\nЛиния №%u\n',i);
    printresult(d(:,i),0,'d','d','пикс');
    printresult(d_prime(:,i),0,'d','d′','пикс');
    printresult(delta_nu(:,i),d_delta_nu(:,i),'pm','dν','м^−1');
    printresult(mean_delta_nu(i),d_mean_delta_nu(i),'pm','mean_dν','м^−1');
    printresult(A(i),d_A(i),'pm','A','Дж','auto');
    printresult(g(i),d_g(i),'pm','g(I)');
    printresult(mu_I_div_mu_nucl(i),d_mu_I_div_mu_nucl(i),'pm','μ_I/μ_яд');
end