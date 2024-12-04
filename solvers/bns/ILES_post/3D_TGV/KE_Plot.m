close all
clear
format long

% Import Reference 1 Solution Data
t_ref  = readtable("DATA.csv",Range="A3:A2002",ReadVariableNames=false);
u_ref  = readtable("DATA.csv",Range="B3:B2002",ReadVariableNames=false);
u_ref  = table2array(u_ref);
t_ref  = table2array(t_ref);

% Import BNS Solution Data
t_soln = readtable("DATA.csv",Range="G3:G2002",ReadVariableNames=false);
u_soln = readtable("DATA.csv",Range="H3:H2002",ReadVariableNames=false);
u_soln  = table2array(u_soln);
t_soln  = table2array(t_soln);

% Import BNS Solution Data
t_soln1 = readtable("DATA.csv",Range="K3:K2002",ReadVariableNames=false);
u_soln1 = readtable("DATA.csv",Range="L3:L2002",ReadVariableNames=false);
u_soln1  = table2array(u_soln1);
t_soln1  = table2array(t_soln1);
 
% Import BNS Solution Data
t_soln2 = readtable("DATA.csv",Range="O3:O2002",ReadVariableNames=false);
u_soln2 = readtable("DATA.csv",Range="P3:P2002",ReadVariableNames=false);
u_soln2  = table2array(u_soln2);
t_soln2  = table2array(t_soln2);

% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="S3:S2002",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="T3:T2002",ReadVariableNames=false);
u_soln3  = table2array(u_soln3);
t_soln3  = table2array(t_soln3);

% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AE3:AE2002",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AF3:AF2002",ReadVariableNames=false);
u_soln4  = table2array(u_soln4);
t_soln4  = table2array(t_soln4);

% Import BNS Solution Data
t_soln5 = readtable("DATA.csv",Range="AI3:AI2002",ReadVariableNames=false);
u_soln5 = readtable("DATA.csv",Range="AJ3:AJ2002",ReadVariableNames=false);
u_soln5  = table2array(u_soln5);
t_soln5  = table2array(t_soln5);

% Plot
hfig = figure;  % save the figure handle in a variable

plot(t_ref,u_ref,'k-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
hold on;
% plot(t_soln,u_soln,'r-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln1,u_soln1,'g-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
% plot(t_soln2,u_soln2,'b-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln3,u_soln3,'m-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
% plot(t_soln4,u_soln4,'b-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln5,u_soln5,'b-','LineWidth',2.3,'DisplayName','$\Omega(t)$');


legend('DNS','BNS-8x4','BNS-16x4','BNS-32x4')


xlabel('t*')
ylabel('$E_{k}$')
fname = 'Uvelocity';

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
hold off


