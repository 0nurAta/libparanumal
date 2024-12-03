clear 
clc
close all

% Read Reference Data
k_ref  = readtable("DATA.csv",Range="E3:E110",ReadVariableNames=false);
Ek_ref  = readtable("DATA.csv",Range="F3:F110",ReadVariableNames=false);
k_ref  = table2array(k_ref(:,1));
Ek_ref  = table2array(Ek_ref(:,1));

% Read Result Data
[k1,E_K1] = EnergySpec3D("fields_32x32_4.txt");
[k2,E_K2] = EnergySpec3D("fields_24x24_4.txt");
[k3,E_K3] = EnergySpec3D("fields_24x24_41.txt");
% [k4,E_K4] = EnergySpec3D("fields_16x16_3.txt");
% [k5,E_K5] = EnergySpec3D("fields_8x8_3.txt");

% Plot
hfig = figure;  % save the figure handle in a variable
loglog(k1, E_K1,'r-', 'LineWidth', 1.5);
hold on;
loglog(k2, E_K2,'b-', 'LineWidth', 1.5);
loglog(k3, E_K3, 'LineWidth', 1.5);
% loglog(k4, E_K4, 'LineWidth', 1.5);
% loglog(k5, E_K5, 'LineWidth', 1.5);
loglog(k_ref, Ek_ref, 'LineWidth', 1.5);
xlabel('Wavenumber k');
ylabel('Energy Spectrum E(k)');
title('Energy Spectrum');
grid on;
fname = 'Energy_Spectrum';
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpng','-painters')