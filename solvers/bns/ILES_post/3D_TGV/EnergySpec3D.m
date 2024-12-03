function [k,E_K] = EnergySpec3D(filename)

% Solution Information
L  = 2*pi;

% % Read Reference Data
% k_ref  = readtable("DATA.csv",Range="E3:E110",ReadVariableNames=false);
% Ek_ref  = readtable("DATA.csv",Range="F3:F110",ReadVariableNames=false);
% k_ref  = table2array(k_ref(:,1));
% Ek_ref  = table2array(Ek_ref(:,1));

% Read Solution Data
soln_data = readFields3D(filename);
u_soln = soln_data(:, 5);
v_soln = soln_data(:, 6);
w_soln = soln_data(:, 7);

n = round((size(u_soln,1))^(1/3));
u = zeros(n,n,n);
v = zeros(n,n,n);
w = zeros(n,n,n);

k = 1;
for i = 1:n
    for j = 1:n
        for l = 1:n
    u(i,j,l) = u_soln(k);
    v(i,j,l) = v_soln(k); 
    w(i,j,l) = w_soln(k); 
    k = k+1;
        end
    end
end

% u_k = fft2(u) / (n * n);  % Normalize
% v_k = fft2(v) / (n * n);  % Normalize

% u_k = fft(u,n,1);  
% u_k = fft(u_k,n,2);
% u_k = fft(u_k,n,3) / (n*n*n); % Normalize
% v_k = fft(v,n,1);  
% v_k = fft(v_k,n,2);
% v_k = fft(v_k,n,3) / (n*n*n); % Normalize
% w_k = fft(w,n,1);  
% w_k = fft(w_k,n,2);
% w_k = fft(w_k,n,3) / (n*n*n); % Normalize

u_k = fftn(u) / (n^3);
v_k = fftn(v) / (n^3);
w_k = fftn(w) / (n^3);
% Compute energy spectrum
E_k = 0.5 * (u_k.*conj(u_k) + v_k.*conj(v_k)+ w_k.*conj(w_k));
E_K = zeros(n, 1);  % Correct size for energy bins

% Compute Wave Numbers
k_x =  fftfreq(n,(L/(n-1)));
k_y =  fftfreq(n,(L/(n-1)));
k_z =  fftfreq(n,(L/(n-1)));
% % 
% k_X = zeros(n,n);
% k_Y = zeros(n,n);
% k_Z = zeros(n,n);
% 
% for i = 1:n
% k_X(i,:) = k_x;
% k_Y(i,:) = k_y;
% k_Z(i,:) = k_z;
% end
% 
% k_X2 = zeros(n,n,n);
% k_Y2 = zeros(n,n,n);
% k_Z2 = zeros(n,n,n);
% 
% for i = 1:n
% k_X2(i,:,:) = k_X;
% k_Y2(i,:,:) = k_Y;
% k_Z2(i,:,:) = k_Z;
% end
% Wave number calculation
k_x = (0:n-1) - floor(n/2);
k_x = (2*pi/L) * fftshift(k_x);
k_x = [zeros(1,1) ,k_x];
k_x = k_x(1,1:n);
[k_X, k_Y, k_Z] = ndgrid(k_x, k_x, k_x);
k_mag = sqrt(k_X.^2 + k_Y.^2+ k_Z.^2);
% k_mag=zeros(n,n);
kMag = zeros(n);
 for j = 1:n
    for i = 1:n
        for l = 1:n
        k_id = round(k_mag(i,j,l))+1;
            if k_id <= n  % Prevent overflow
            E_K(k_id) = E_K(k_id) + E_k(i, j,l);  % Accumulate enery
            end
        end
    end
 end
k = ((0:n-1))*(2*pi/L);  % Wavenumber range
% Plot the energy spectrum

% % Plot
% hfig = figure;  % save the figure handle in a variable
% loglog(k, E_K, 'LineWidth', 1.5);
% hold on;
% loglog(k_ref, Ek_ref, 'LineWidth', 1.5);
% xlabel('Wavenumber k');
% ylabel('Energy Spectrum E(k)');
% title('Energy Spectrum');
% grid on;
% fname = 'Energy_Spectrum';
% picturewidth = 20; % set this parameter and keep it forever
% hw_ratio = 0.65; % feel free to play with this ratio
% set(findall(hfig,'-property','FontSize'),'FontSize',15) % adjust fontsize to your document
% set(findall(hfig,'-property','Box'),'Box','off') % optional
% set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(hfig,'Position');
% set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% %print(hfig,fname,'-dpdf','-painters','-fillpage')
% print(hfig,fname,'-dpng','-painters')

end