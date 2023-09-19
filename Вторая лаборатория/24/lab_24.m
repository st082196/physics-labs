load('data_24');

% λ, нм
wl_peaks = wl(N_peaks);
wl_g = wl(N_g);
% E, эВ
E = h.*c./(n_0.*wl.*1e-9.*e);
E_peaks = E(N_peaks);
d_E_peaks = abs(E_peaks.*d_wl_peaks./wl_peaks);
E_g = E(N_g);
d_E_g = abs(E_g.*d_wl_g./wl_g);
% Полуширина, мэВ
FWHM = zeros(size(N_peaks));
for i = 1:numel(N_peaks)
    peak = peak_bounds(i,1):peak_bounds(i,2);
    FWHM(i) = peakwidth(E(peak),I(peak)).*1e3;
end
% L, нм
L_nm = h./2.*sqrt((m_e.^-1 + m_h.^-1)./(2.*(E_peaks - E_g).*e)).*1e9;
d_L_nm = abs(L_nm./(2.*(E_peaks - E_g)).*sqrt(d_E_peaks.^2 + d_E_g.^2));
% L, мс
L_ml = L_nm./(a_0./2);
d_L_ml = d_L_nm./(a_0./2);
%% Вывод результатов
printresult(wl_g,d_wl_g,'()','λ_g','нм');
printresult(E_g,d_E_g,'()','E_g','эВ');
printresult(wl_peaks,d_wl_peaks,'()','λ','нм');
printresult(E_peaks,d_E_peaks,'()','E','эВ');
printresult(FWHM,3,'s','Полуширина','мэВ');
printresult(L_nm,d_L_nm,'()','L','нм');
printresult(L_ml,d_L_ml,'()','L','мс');
%% I(λ)
figure('Name','I(λ)','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:size(I,2)
    plotobject = plot(wl(:,i),I(:,i),'Color',color{i});
    plotobject.DataTipTemplate.DataTipRows(3).Label = 'Index';
    plotobject.DataTipTemplate.DataTipRows(3).Value = 1:numel(plotobject.XData);
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itλ\rm, нм','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
xlim([min(wl) max(wl)]);
pbaspect([2 1 1]);
figresize('width',17);
%% I(E)
figure('Name','I(E)','Units','centimeters');
hold on;

color = {'b','r','g'};
for i = 1:size(I,2)
    plot(E(:,i),I(:,i),'Color',color{i});
end

set(gca,'FontName','Calibri','FontSize',14,'Box','on','LooseInset',get(gca,'TightInset'));
xlabel('\itE\rm, эВ','FontName','Liberation Serif');
ylabel('\itI','FontName','Liberation Serif');
xlim([min(E) max(E)]);
pbaspect([2 1 1]);
figresize('width',17);