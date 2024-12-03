close all
clear
format long

% 
T = 0:0.01:20;
T1 = 0:0.01:20;
% Import Reference 1 Solution Data
u_ref1  = readtable("Ref_Solution.csv",Range="F3:F2002",ReadVariableNames=false);
t_ref1  = readtable("Ref_Solution.csv",Range="D3:D2002",ReadVariableNames=false);
u_ref1  = table2array(u_ref1);
t_ref1  = table2array(t_ref1);

x = lsqr(u_ref1,t_ref1);
% Import BNS Solution Data
u_soln = readtable("Ref_Solution.csv",Range="AK3:AK202",ReadVariableNames=false);
t_soln = readtable("Ref_Solution.csv",Range="AJ3:AJ202",ReadVariableNames=false);
u_soln = table2array(u_soln);
t_soln = table2array(t_soln);

u_soln = spline(t_soln,u_soln, T1);
u_soln = -gradient(u_soln, T);

% Import BNS Solution Data
u_soln1 = readtable("Ref_Solution.csv",Range="AQ3:AQ2002",ReadVariableNames=false);
t_soln1 = readtable("Ref_Solution.csv",Range="AP3:AP2002",ReadVariableNames=false);
u_soln1 = table2array(u_soln1);
t_soln1 = table2array(t_soln1);

u_soln1 = spline(t_soln1,u_soln1, T);
u_soln1 = -gradient(u_soln1, T);

%Import BNS Solution Data
u_soln2 = readtable("Ref_Solution.csv",Range="AU3:AU2002",ReadVariableNames=false);
t_soln2 = readtable("Ref_Solution.csv",Range="AT3:AT2002",ReadVariableNames=false);
u_soln2 = table2array(u_soln2);
t_soln2 = table2array(t_soln2);

u_soln2 = spline(t_soln2,u_soln2, T);
u_soln2 = -gradient(u_soln2, T);

%Import BNS Solution Data
u_soln3 = readtable("Ref_Solution.csv",Range="BC3:BC2002",ReadVariableNames=false);
t_soln3 = readtable("Ref_Solution.csv",Range="BB3:BB2002",ReadVariableNames=false);
u_soln3 = table2array(u_soln3);
t_soln3 = table2array(t_soln3);

u_soln3 = spline(t_soln3,u_soln3, T);
u_soln3 = -gradient(u_soln3, T);

% Plot
hfig = figure;  % save the figure handle in a variable
figure
plot(t_ref1,u_ref1,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
hold on;
plot(T,u_soln,'r-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(T,u_soln1,'g-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(T,u_soln2,'b-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(T,u_soln3,'m-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% plot(K,E3,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');

legend('DNS','BNS-8x8','BNS-16x16','BNS-24x24')


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



% % 
% pp = spline(t_soln3,u_soln3, 0:0.3:20);
% % piece-wise polynomial
% %[breaks,coefs,l,k,d] = unmkpp(pp);
% % get its derivative
% %pp_der= diff(pp, 0:0.5:20);
% %pp_der = mkpp(breaks,repmat(k-1:-1:1,d*l,1).*coefs(:,1:k-1),d);
% pp_der = -gradient(pp, 0:0.3:20);
% %a = -ppval(pp_der,0:0.5:20);