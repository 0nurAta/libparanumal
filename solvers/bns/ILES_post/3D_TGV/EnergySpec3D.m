function [k,E_K] = EnergySpec3D(filename)
  format long
  
% Solution Information
L  = 2*pi;

% % Read Reference Data
% k_ref  = readtable("DATA.csv",Range="E3:E110",ReadVariableNames=false);
% Ek_ref  = readtable("DATA.csv",Range="F3:F110",ReadVariableNames=false);
% k_ref  = table2array(k_ref(:,1));
% Ek_ref  = table2array(Ek_ref(:,1));

% Read Solution Data
soln_data = readFields3D(filename);
%soln_data = readmatrix(filename);
u_soln = soln_data(:, 4);
v_soln = soln_data(:, 5);
w_soln = soln_data(:, 6);

n = round((size(u_soln,1))^(1/3));
u = zeros(n,n,n);
v = zeros(n,n,n);
w = zeros(n,n,n);

% k = 1;
% for i = 1:n
%     for j = 1:n
%         for l = 1:n
%     u(i,j,l) = u_soln(k);
%     v(i,j,l) = v_soln(k); 
%     w(i,j,l) = w_soln(k); 
%     k = k+1;
%         end
%     end
% end
u = reshape(u_soln, n, n, n);
v = reshape(v_soln, n, n, n);
w = reshape(w_soln, n, n, n);
% u_k = fft2(u) / (n * n);  % Normalize
% v_k = fft2(v) / (n * n);  % Normalize

u_k = fft(u,n,1);  
u_k = fft(u_k,n,2);
u_k = fft(u_k,n,3) / (n*n*n); % Normalize
v_k = fft(v,n,1);  
v_k = fft(v_k,n,2);
v_k = fft(v_k,n,3) / (n*n*n); % Normalize
w_k = fft(w,n,1);  
w_k = fft(w_k,n,2);
w_k = fft(w_k,n,3) / (n*n*n); % Normalize

% u_k = (fftn(u) / (n^3));
% v_k = fftn(v) / (n^3);
% w_k = fftn(w) / (n^3);


% Compute Wave Numbers
k_x =  fftfreq(n,(L/(n-1))).*(2*pi);
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
% k_x = (0:n-1) - floor(n/2);
% k_x = (2*pi/L) * fftshift(k_x);
% k_x = [zeros(1,1) ,k_x];
% k_x = k_x(1,1:n);


[k_X, k_Y, k_Z] = ndgrid(k_x, k_x, k_x);
k_mag = sqrt(k_X.^2 + k_Y.^2+ k_Z.^2);
dk = k_x(2)-k_x(1);
% Compute energy spectrum
% E_k = 0.5 * (u_k.*conj(u_k) + v_k.*conj(v_k)+ w_k.*conj(w_k))/dk;
% E_ku = u_k.*conj(u_k);
% E_kv = v_k.*conj(v_k);
% E_kw = w_k.*conj(w_k);
E_ku = fftshift(u_k.*conj(u_k));
E_kv = fftshift(v_k.*conj(v_k));
E_kw = fftshift(w_k.*conj(w_k));

centerx = floor(n/2.);
centery = floor(n/2.);
centerz = floor(n/2.);
rad_sphere = ceil(sqrt((3*(n*n)))/2.+1);

epsilon = 1e-50;

E_Ku = zeros(rad_sphere, 1)+epsilon;  % Correct size for energy bins
E_Kv = zeros(rad_sphere, 1)+epsilon;  % Correct size for energy bins
E_Kw = zeros(rad_sphere, 1)+epsilon;  % Correct size for energy bins




 for j = 1:rad_sphere
    for i = 1:rad_sphere
        for l = 1:rad_sphere
        % k_id = round(k_mag(i,j,l))+1;
        k_id =  round(sqrt((i-centerx)^2+(j-centery)^2+(l-centerz)^2))+1;
             if k_id <= rad_sphere %&& k_id > 0.1  % Prevent overflow
            % E_K(k_id) = E_K(k_id) + E_k(i, j,l);  % Accumulate enery
            E_Ku(k_id) = E_Ku(k_id) + E_ku(i, j,l);  % Accumulate enery
            E_Kv(k_id) = E_Kv(k_id) + E_kv(i, j,l);  % Accumulate enery
            E_Kw(k_id) = E_Kw(k_id) + E_kw(i, j,l);  % Accumulate enery
            end
        end
    end
 end

E_K = 0.5 * (E_Ku+E_Kv+E_Kw)/dk;
k = (((0:rad_sphere-1))*(2*pi)/L);  % Wavenumber range

E_K = E_K(3:rad_sphere-4,1);
k = k(1,3:rad_sphere-4);
end