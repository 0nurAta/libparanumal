close all
clear
format long

% 
T = 0:0.01:20;
T1 = 0:0.01:20;
% Import Reference 1 Solution Data
t_ref  = readtable("DATA.csv",Range="A3:A2002",ReadVariableNames=false);
u_ref  = readtable("DATA.csv",Range="D3:D2002",ReadVariableNames=false);
u_ref  = table2array(u_ref);
t_ref  = table2array(t_ref);

%
% Import BNS Solution Data
t_soln2 = readtable("DATA.csv",Range="O3:O2002",ReadVariableNames=false);
u_soln2 = readtable("DATA.csv",Range="R3:R2002",ReadVariableNames=false);
u_soln2  = -table2array(u_soln2(1:1:2000,1));
t_soln2  = table2array(t_soln2(1:1:2000,1));

    % u_soln2 = spline(t_soln2,u_soln2, T);
% u_soln2 = -gradient(u_soln2, T);

% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="W3:W2002",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="Z3:Z2002",ReadVariableNames=false);
u_soln3  = -table2array(u_soln3(1:1:2000,1));
t_soln3  = table2array(t_soln3(1:1:2000,1));
 
% Sample energy data (replace this with your actual data)
%time = 0:0.01:20; % Time array

% u_soln3 = spline(t_soln3,u_soln3, T);
% u_soln3 = -gradient(u_soln3, T);

% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AE3:AE2002",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AF3:AF2002",ReadVariableNames=false);
u_soln4  = table2array(u_soln4(1:1:2000,1));
t_soln4  = table2array(t_soln4(1:1:2000,1));

u_soln4 = spline(t_soln4,u_soln4, T);
u_soln4 = -gradient(u_soln4, T);

% Import BNS Solution Data
t_soln5 = readtable("DATA.csv",Range="AM3:AM2002",ReadVariableNames=false);
u_soln5 = readtable("DATA.csv",Range="AN3:AN2002",ReadVariableNames=false);
u_soln5  = table2array(u_soln5(1:1:2000,1));
t_soln5  = table2array(t_soln5(1:1:2000,1));

u_soln5 = spline(t_soln5,u_soln5, T);
u_soln5 = -gradient(u_soln5, T);

% Plot
hfig = figure;  % save the figure handle in a variable 
plot(t_soln2,u_soln2,'b-','LineWidth',1.0);hold on;
plot(t_soln3,u_soln3,'c--','LineWidth',1.0);
plot(T,u_soln4,'m-.','LineWidth',1.0);
plot(T,u_soln5,'g:','LineWidth',1.5);
plot(t_ref,u_ref,'k-','LineWidth',2.0);
xlim([0, 20]);
ylim([0.5e-3, 14.5e-3]);
% xlim([5, 15]);
% ylim([7e-3, 14e-3]);
% Set Legends
string1 = '$$BNS-8\times3$$';
string2 = '$$BNS-16\times3$$';
string3 = '$$BNS-24\times3$$';
string4 = '$$BNS-32\times3$$';
string5 = '$$BNS-48\times3$$';
string6 = '$$Spectral-512$$';

leg1 = legend(string2,string3,string4,string5,string6);



picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',22) % adjust fontsize to your document
fs=18;
set(leg1,'interpreter','latex','fontsize',fs,'Location','NorthEast','box','on');
xlabel('$$t$$','interpreter','latex','fontsize',fs+4);
ylabel('$$\epsilon_{1}$$','interpreter','latex','fontsize',fs+4);
fname = 'Epsilon1';
set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hfig,fname,'-dpdf','-painters','-fillpage')
%print(hfig,fname,'-dpng','-painters')
print(hfig,fname,'-deps','-painters')
hold off
