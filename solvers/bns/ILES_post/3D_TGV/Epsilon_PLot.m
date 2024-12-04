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
t_soln = readtable("DATA.csv",Range="G3:G2002",ReadVariableNames=false);
u_soln = readtable("DATA.csv",Range="H3:H2002",ReadVariableNames=false);
u_soln  = table2array(u_soln);
t_soln  = table2array(t_soln);

% u_soln = spline(t_soln,u_soln, T1);
% u_soln = -gradient(u_soln, T);


% Import BNS Solution Data
t_soln1 = readtable("DATA.csv",Range="K3:K2002",ReadVariableNames=false);
u_soln1 = readtable("DATA.csv",Range="L3:L2002",ReadVariableNames=false);
u_soln1  = table2array(u_soln1);
t_soln1  = table2array(t_soln1);
 
% Sample energy data (replace this with your actual data)
%time = 0:0.01:20; % Time array

u_soln1 = spline(t_soln1,u_soln1, T);
u_soln1 = -gradient(u_soln1, T);

% Define moving average window size
%window_size = 250;

% Moving average smoothing
%u_soln1 = movmean(u_soln1, window_size);

% Import BNS Solution Data
t_soln2 = readtable("DATA.csv",Range="O3:O2002",ReadVariableNames=false);
u_soln2 = readtable("DATA.csv",Range="P3:P2002",ReadVariableNames=false);
u_soln2  = table2array(u_soln2);
t_soln2  = table2array(t_soln2);

u_soln2 = spline(t_soln2,u_soln2, T);
u_soln2 = -gradient(u_soln2, T);

% Import BNS Solution Data
t_soln3 = readtable("DATA.csv",Range="S3:S2002",ReadVariableNames=false);
u_soln3 = readtable("DATA.csv",Range="T3:T2002",ReadVariableNames=false);
u_soln3  = table2array(u_soln3);
t_soln3  = table2array(t_soln3);

u_soln3 = spline(t_soln3,u_soln3, T);
u_soln3 = -gradient(u_soln3, T);

% Import BNS Solution Data
t_soln4 = readtable("DATA.csv",Range="AE3:AE2002",ReadVariableNames=false);
u_soln4 = readtable("DATA.csv",Range="AF3:AF2002",ReadVariableNames=false);
u_soln4  = table2array(u_soln4);
t_soln4  = table2array(t_soln4);

u_soln4 = spline(t_soln4,u_soln4, T);
u_soln4 = -gradient(u_soln4, T);

% Import BNS Solution Data
t_soln5 = readtable("DATA.csv",Range="AI3:AI2002",ReadVariableNames=false);
u_soln5 = readtable("DATA.csv",Range="AJ3:AJ2002",ReadVariableNames=false);
u_soln5  = table2array(u_soln5);
t_soln5  = table2array(t_soln5);

u_soln5 = spline(t_soln5,u_soln5, T);
u_soln5 = -gradient(u_soln5, T);

% Plot
hfig = figure;  % save the figure handle in a variable
plot(t_ref,u_ref,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
hold on;
plot(T,u_soln1,'r-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(T,u_soln3,'g-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
plot(T,u_soln5,'b-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
%plot(T,u_soln3,'m-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% plot(K,E3,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');

legend('DNS','BNS-8x3','BNS-16x4','BNS-32x4')


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