close all
clear
format long
% Import Reference Data
u_ref = readtable("x_-0.48.csv",Range="B2:B91",ReadVariableNames=false);
y_ref = readtable("x_-0.48.csv",Range="C2:C91",ReadVariableNames=false);
u_ref = table2array(u_ref);
y_ref = table2array(y_ref);

% Import Solution Data
u_soln = readtable("n=2.csv",Range="A2:A1002",ReadVariableNames=false);
y_soln = readtable("n=2.csv",Range="B2:B1002",ReadVariableNames=false);
u_soln = table2array(u_soln);
y_soln = table2array(y_soln);

% Import Solution Data
u_soln1 = readtable("n=3.csv",Range="A2:A1002",ReadVariableNames=false);
y_soln1 = readtable("n=3.csv",Range="B2:B1002",ReadVariableNames=false);
u_soln1 = table2array(u_soln1);
y_soln1 = table2array(y_soln1);

% Plot
hfig = figure;  % save the figure handle in a variable

plot(u_ref,y_ref,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
hold on;
plot(u_soln,y_soln,'r-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(u_soln1,y_soln1,'b-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% plot(K,E3,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');

legend('DNS','BNS-p=2','BNS-p=3')


xlabel('x')
ylabel('u')
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
