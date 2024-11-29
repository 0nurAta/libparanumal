close all
clear
format long

% Import Reference 1 Solution Data
u_ref1  = readtable("Ref_Solution.csv",Range="A3:A59",ReadVariableNames=false);
t_ref1  = readtable("Ref_Solution.csv",Range="B3:B59",ReadVariableNames=false);
u_ref1  = table2array(u_ref1);
t_ref1  = table2array(t_ref1);

% Import BNS Solution Data
u_soln = readtable("Ref_Solution.csv",Range="AK3:AK202",ReadVariableNames=false);
t_soln = readtable("Ref_Solution.csv",Range="AJ3:AJ202",ReadVariableNames=false);
u_soln = table2array(u_soln);
t_soln = table2array(t_soln);

% Import BNS Solution Data
u_soln1 = readtable("Ref_Solution.csv",Range="AQ3:AQ2002",ReadVariableNames=false);
t_soln1 = readtable("Ref_Solution.csv",Range="AP3:AP2002",ReadVariableNames=false);
u_soln1 = table2array(u_soln1);
t_soln1 = table2array(t_soln1);
% 
%Import BNS Solution Data
u_soln2 = readtable("Ref_Solution.csv",Range="AU3:AU2002",ReadVariableNames=false);
t_soln2 = readtable("Ref_Solution.csv",Range="AT3:AT2002",ReadVariableNames=false);
u_soln2 = table2array(u_soln2);
t_soln2 = table2array(t_soln2);

%Import BNS Solution Data
u_soln3 = readtable("Ref_Solution.csv",Range="BC3:BC2002",ReadVariableNames=false);
t_soln3 = readtable("Ref_Solution.csv",Range="BB3:BB2002",ReadVariableNames=false);
u_soln3 = table2array(u_soln3);
t_soln3 = table2array(t_soln3);

% Plot
hfig = figure;  % save the figure handle in a variable

plot(t_ref1,u_ref1,'k-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
hold on;
plot(t_soln,u_soln,'r-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln1,u_soln1,'g-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln2,u_soln2,'b-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
plot(t_soln3,u_soln3,'m-','LineWidth',2.3,'DisplayName','$\Omega(t)$');
% plot(K,E3,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');

legend('DNS','BNS-8x8','BNS-16x16','BNS-24x24','BNS-32x32')


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


