clear 
clc
close all

% Read Reference Data
k_ref  = readtable("DATA.csv",Range="E3:E110",ReadVariableNames=false);
Ek_ref  = readtable("DATA.csv",Range="F3:F110",ReadVariableNames=false);
k_ref  = table2array(k_ref(:,1));
Ek_ref  = table2array(Ek_ref(:,1));

% Read Reference Data
data1 = readmatrix("Ek_Spectrum_16x3.dat");
data2 = readmatrix("Ek_Spectrum.dat");
data3 = readmatrix("Ek_Spectrum_32x3.dat");
data4 = readmatrix("Ek_Spectrum_48x3.dat");

% Attain Data
k1   = data1(3:(size(data1,1)-5),1);
E_K1 = data1(3:(size(data1,1)-5),2);
k2   = data2(3:(size(data2,1)-5),1);
E_K2 = data2(3:(size(data2,1)-5),2);
k3   = data3(3:(size(data3,1)-6),1);
E_K3 = data3(3:(size(data3,1)-6),2);
k4   = data4(3:(size(data4,1)-10),1);
E_K4 = data4(3:(size(data4,1)-10),2);
% k5   = data5(3:(size(data5,1)-1),1);
% E_K5 = data5(3:(size(data5,1)-1),2);

% Read Result Data
% [k1,E_K1] = EnergySpec3D("fields_32x32_3.txt");
% [k2,E_K2] = EnergySpec3D("fields_8x8_3.txt");
% [k3,E_K3] = EnergySpec3D("fields_24x24_41.txt");
% [k4,E_K4] = EnergySpec3D("fields_16x16_3.txt");
% [k5,E_K5] = EnergySpec3D("fields_8x8_3.txt");

% Plot Results
hfig = figure;  % save the figure handle in a variable

loglog(k1, E_K1,'b-', 'MarkerSize',2, 'LineWidth', 1.0); hold on;
loglog(k2, E_K2,'c--', 'LineWidth', 1.0);
loglog(k3, E_K3,'m-.', 'LineWidth', 1.0);
loglog(k4, E_K4,'g:', 'LineWidth', 1.5);
loglog(k_ref, Ek_ref,'k-' ,'LineWidth', 2.00); 
xlim([10e-1, 5e+2]);
ylim([10e-11, 10e-1]);
xline(24,'-k');
xline(36,'--k');
xline(48,'-.k');
xline(72,':k');


% Set Legends
% string1 = '$$BNS-8\times3$$';
string2 = '$$BNS-16\times3$$';
string3 = '$$BNS-24\times3$$';
string4 = '$$BNS-32\times3$$';
string5 = '$$BNS-48\times3$$';
string6 = '$$Spectral$$';
leg = legend( string2,string3,string4,string5,string6);
%xlabel('Wavenumber k');
%ylabel('Energy Spectrum E(k)');
fs=18;
set(leg,'interpreter','latex','fontsize',fs+2,'Location','SouthWest','box','off');
xlabel('$$k$$','interpreter','latex','fontsize',fs+2);
ylabel('$$E(k)$$','interpreter','latex','fontsize',fs+2);

%grid on;
fname = 'Energy_Spectrum';
picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',22) % adjust fontsize to your document
set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
set(gca,'YMinorTick','off')
print(hfig,fname,'-deps','-painters')
%print(hfig,fname,'-dpng','-painters')