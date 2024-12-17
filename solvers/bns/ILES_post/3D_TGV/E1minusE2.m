close all
clear
format long

% Import BNS Solution Data
t_soln = readtable("DATA.csv",Range="AU3:AU2001",ReadVariableNames=false);
u_soln = readtable("DATA.csv",Range="AX3:AX2001",ReadVariableNames=false);
t_soln  = table2array(t_soln);
u_soln  = -table2array(u_soln);
% Import BNS Solution Data
t_soln1 = readtable("DATA.csv",Range="AU3:AU2001",ReadVariableNames=false);
u_soln1 = readtable("DATA.csv",Range="AW3:AW2001",ReadVariableNames=false);
t_soln1  = table2array(t_soln1);
u_soln1  = table2array(u_soln1);

u_soln1 = 2*0.000625*u_soln1;

t_soln2 = t_soln1;
u_soln2 = u_soln-u_soln1;

% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="AQ3:AQ2001",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="AT3:AT2001",ReadVariableNames=false);
t_soln3  = table2array(t_soln3);
u_soln3  = -table2array(u_soln3);
% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AQ3:AQ2001",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AS3:AS2001",ReadVariableNames=false);
t_soln4  = table2array(t_soln4);
u_soln4  = table2array(u_soln4);

u_soln4 = 2*0.000625*u_soln4;

t_soln5 = t_soln4;
u_soln5 = u_soln3-u_soln4;

% Import BNS Solution Data
t_soln6 = readtable("DATA.csv",Range="AI3:AI2001",ReadVariableNames=false);
u_soln6 = readtable("DATA.csv",Range="AL3:AL2001",ReadVariableNames=false);
t_soln6  = table2array(t_soln6);
u_soln6  = -table2array(u_soln6);
% Import BNS Solution Data
t_soln7 = readtable("DATA.csv",Range="AI3:AI2001",ReadVariableNames=false);
u_soln7 = readtable("DATA.csv",Range="AK3:AK2001",ReadVariableNames=false);
t_soln7  = table2array(t_soln7);
u_soln7  = table2array(u_soln7);

u_soln7 = 2*0.000625*u_soln7;

t_soln8 = t_soln7;
u_soln8 = u_soln6-u_soln7;

% Import BNS Solution Data
t_soln9 = readtable("DATA.csv",Range="AE3:AE2001",ReadVariableNames=false);
u_soln9 = readtable("DATA.csv",Range="AH3:AH2001",ReadVariableNames=false);
t_soln9  = table2array(t_soln9);
u_soln9  = -table2array(u_soln9);
% Import BNS Solution Data
t_soln10 = readtable("DATA.csv",Range="AE3:AE2001",ReadVariableNames=false);
u_soln10 = readtable("DATA.csv",Range="AG3:AG2001",ReadVariableNames=false);
t_soln10  = table2array(t_soln10);
u_soln10  = table2array(u_soln10);

u_soln10 = 2*0.000625*u_soln10;

t_soln11 = t_soln10;
u_soln11 = u_soln9-u_soln10;

% Plot
hfig = figure;  % save the figure handle in a variable 
plot(t_soln11,u_soln11,'b-','LineWidth',1.0);hold on;
plot(t_soln8,u_soln8,'c--','LineWidth',1.0); hold on;
plot(t_soln5,u_soln5,'m-.','LineWidth',1.0);
plot(t_soln2,u_soln2,'g:','LineWidth',1.5);
%plot(t_ref,u_ref,'k-','LineWidth',2.0);
% xlim([-0, 20]);
ylim([-0.0005, 0.002]);
% xlim([7, 12]);
% ylim([8e-3, 14e-3]);
% Set Legends
string1 = '$$BNS-32\times3$$';
string2 = '$$BNS-32\times4$$';
string3 = '$$BNS-32\times5$$';
string4 = '$$BNS-32\times6$$';
% 
leg1 = legend(string1, string2, string3, string4);



picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',22) % adjust fontsize to your document
fs=18;
set(leg1,'interpreter','latex','fontsize',fs+2,'Location','NorthEast','box','on');
xlabel('$$t$$','interpreter','latex','fontsize',fs+4);
ylabel('$$\epsilon_{num}$$','interpreter','latex','fontsize',fs+4);
fname = 'Enstrophy';
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


