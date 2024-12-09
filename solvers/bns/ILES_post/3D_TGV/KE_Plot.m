close all
clear
format long

% Import Reference 1 Solution Data
t_ref  = readtable("DATA.csv",Range="A3:A2002",ReadVariableNames=false);
u_ref  = readtable("DATA.csv",Range="B3:B2002",ReadVariableNames=false);
u_ref  = table2array(u_ref);
t_ref  = table2array(t_ref);

% Import BNS Solution Data
% t_soln1 = readtable("DATA.csv",Range="G3:G2002",ReadVariableNames=false);
% u_soln1 = readtable("DATA.csv",Range="H3:H2002",ReadVariableNames=false);
% u_soln1  = table2array(u_soln1);
% t_soln1  = table2array(t_soln1);
 
% Import BNS Solution Data
t_soln2 = readtable("DATA.csv",Range="O3:O2002",ReadVariableNames=false);
u_soln2 = readtable("DATA.csv",Range="P3:P2002",ReadVariableNames=false);
u_soln2  = table2array(u_soln2);
t_soln2  = table2array(t_soln2);

% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="W3:W2002",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="X3:X2002",ReadVariableNames=false);
u_soln3  = table2array(u_soln3);
t_soln3  = table2array(t_soln3);

% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AE3:AE2002",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AF3:AF2002",ReadVariableNames=false);
u_soln4  = table2array(u_soln4);
t_soln4  = table2array(t_soln4);

% Import BNS Solution Data
t_soln5 = readtable("DATA.csv",Range="AM3:AM2002",ReadVariableNames=false);
u_soln5 = readtable("DATA.csv",Range="AN3:AN2002",ReadVariableNames=false);
u_soln5  = table2array(u_soln5);
t_soln5  = table2array(t_soln5);

% Plot
hfig = figure;  % save the figure handle in a variable 
plot(t_soln2,u_soln2,'b-','LineWidth',1.0);hold on;
plot(t_soln3,u_soln3,'c--','LineWidth',1.0);
plot(t_soln4,u_soln4,'m-.','LineWidth',1.0);
plot(t_soln5,u_soln5,'g:','LineWidth',1.5);
plot(t_ref,u_ref,'k-','LineWidth',2.0);
% xlim([10,20]);
% Set Legends
string1 = '$$Spectral-512$$';
string2 = '$$BNS-8\times3$$';
string3 = '$$BNS-16\times3$$';
string4 = '$$BNS-24\times3$$';
string5 = '$$BNS-32\times3$$';
string6 = '$$BNS-48\times3$$';
leg1 = legend( string3,string4,string5,string6,string1);

fs=18;
set(leg1,'interpreter','latex','fontsize',fs+2,'Location','NorthEast','box','on');
xlabel('$$t$$','interpreter','latex','fontsize',fs+4);
ylabel('$$E_{k}$$','interpreter','latex','fontsize',fs+4);
fname = 'Kinetic_Energy';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',22) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','on') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-deps','-painters')
hold off


