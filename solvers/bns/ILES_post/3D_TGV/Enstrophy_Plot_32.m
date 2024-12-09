close all
clear
format long

% Import Reference 1 Solution Data
t_ref  = readtable("DATA.csv",Range="A3:A2002",ReadVariableNames=false);
u_ref  = readtable("DATA.csv",Range="C3:C2002",ReadVariableNames=false);
u_ref  = table2array(u_ref);
t_ref  = table2array(t_ref);

% % Import BNS Solution Data
% t_soln1 = readtable("DATA.csv",Range="G3:G2002",ReadVariableNames=false);
% u_soln1 = readtable("DATA.csv",Range="I3:I2002",ReadVariableNames=false);
% u_soln1  = table2array(u_soln1);
% t_soln1  = table2array(t_soln1);
 
% Import BNS Solution Data
t_soln2 = readtable("DATA.csv",Range="AE3:AE2002",ReadVariableNames=false);
u_soln2 = readtable("DATA.csv",Range="AG3:AG2002",ReadVariableNames=false);
t_soln2  = table2array(t_soln2);
u_soln2  = table2array(u_soln2);


% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="AI3:AI2002",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="AK3:AK2002",ReadVariableNames=false);
t_soln3  = table2array(t_soln3);
u_soln3  = table2array(u_soln3);


% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AQ3:AQ2002",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AS3:AS2002",ReadVariableNames=false);
t_soln4  = table2array(t_soln4);
u_soln4  = table2array(u_soln4);


% Import BNS Solution Data
t_soln5 = readtable("DATA.csv",Range="AU3:AU2002",ReadVariableNames=false);
u_soln5 = readtable("DATA.csv",Range="AW3:AW2002",ReadVariableNames=false);
u_soln5  = table2array(u_soln5);
t_soln5  = table2array(t_soln5);

u_ref   = 2*0.000625*u_ref;
u_soln2 = 2*0.000625*u_soln2;
u_soln3 = 2*0.000625*u_soln3;
u_soln4 = 2*0.000625*u_soln4;
u_soln5 = 2*0.000625*u_soln5;

% Plot
hfig = figure;  % save the figure handle in a variable 
plot(t_soln2,u_soln2,'b-','LineWidth',1.0);hold on;
plot(t_soln3,u_soln3,'c--','LineWidth',1.0);
plot(t_soln4,u_soln4,'m-.','LineWidth',1.0);
plot(t_soln5,u_soln5,'g:','LineWidth',1.5);
plot(t_ref,u_ref,'k-','LineWidth',2.0);
% xlim([0, 20]);
% ylim([0.5e-3, 14.5e-3]);
xlim([7, 12]);
ylim([8e-3, 14e-3]);
% Set Legends
string1 = '$$BNS-32\times3$$';
string2 = '$$BNS-32\times4$$';
string3 = '$$BNS-32\times5$$';
string4 = '$$BNS-32\times6$$';
string5 = '$$Spectral-512$$';

leg1 = legend(string1,string2,string3,string4,string5);



picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.75; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',22) % adjust fontsize to your document
fs=18;
set(leg1,'interpreter','latex','fontsize',fs+2,'Location','NorthEast','box','on');
xlabel('$$t$$','interpreter','latex','fontsize',fs+4);
ylabel('$$\epsilon_{2}$$','interpreter','latex','fontsize',fs+4);
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


