close all
clear
format long
% Flow variables
L  = pi;    %Domain Length
nu = 0.000625;  %Viscosity
t  = 10;    % Time

% Calculate the exact velocity
% For x velocity
x = -0:0.001:2*L;
y = 3.14159265359*ones(size(x));

% For y velocity
% y = -0:0.001:2*L;
% x = 3.14159265359*ones(size(y));

uexact = sin(x).*cos(y)*exp(-2*nu*t);
vexact = -cos(x).*sin(y)*exp(-2*nu*t);
% Calculate the volume averaged KE

KE = 0.5*9.8696*exp(-2*nu*t)/(6.28318530718)^2

% Calculate Vorticity

w = -cos(x).*cos(y)*exp(-2*nu*t)+sin(x).*sin(y)*exp(-2*nu*t);

% Enstrophy

w = -cos(x).*cos(y)*exp(-2*nu*t)+sin(x).*sin(y)*exp(-2*nu*t);
w2= 0.5*39.478*exp(-2*nu*t)/(6.28318530718)^2
w3D= 0.5*(31.006)*exp(-2*nu*t)/(6.28318530718)^3;
% Import Solution Data
v_soln = readtable("n=8_8x8.csv",Range="C2:C1000",ReadVariableNames=false);
y_soln = readtable("n=8_8x8.csv",Range="H2:H1000",ReadVariableNames=false);
v_soln = table2array(v_soln);
y_soln = table2array(y_soln);
% u_soln = 20.*u_soln;

% Relative Error of Maximum Velocity

Epsilon = abs(max(vexact)-max(v_soln))/max(vexact)


% Plot
hfig = figure;  % save the figure handle in a variable
plot(x,uexact,'k--','LineWidth',2.5,'DisplayName','$\Omega(t)$');
hold on;
plot(y_soln,v_soln,'r:','LineWidth',2.5,'DisplayName','$\Omega(t)$');

% plot(K,E3,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');

legend('Analytical','BNS','Location','northwest')


xlabel('x')
ylabel('u')
fname = 'Vvelocity_TGV2D';

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
set(findall(hfig,'-property','FontSize'),'FontSize',24) % adjust fontsize to your document

set(findall(hfig,'-property','Box'),'Box','off') % optional
set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(hfig,'Position');
set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
%print(hfig,fname,'-dpdf','-painters','-fillpage')
print(hfig,fname,'-dpdf','-painters')
