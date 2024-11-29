% Solution Information
L  = 2*pi;

% Import Solution Data
u_soln = readtable("deneme.csv",Range="A2:A25922",ReadVariableNames=false);
v_soln = readtable("deneme.csv",Range="B2:B25922",ReadVariableNames=false);
u_soln = table2array(u_soln(:,1));
v_soln = table2array(v_soln(:,1));

n = round(sqrt(size(u_soln,1)));
u = zeros(n,n);
v = zeros(n,n);
k = 1;
for i = 1:n
    for j = 1:n
    u(i,j) = u_soln(k);
    v(i,j) = v_soln(k); 
    k = k+1;
    end
end

% u_k = fft2(u) / (n * n);  % Normalize
% v_k = fft2(v) / (n * n);  % Normalize

u_k = fft(u,n,1);  % Normalize
u_k = fft(u_k,n,2) / (n * n);
v_k = fft(v,n,1);  % Normalize
v_k = fft(v_k,n,2) / (n * n);
% u_k = fftshift(u_k);
% v_k = fftshift(v_k);
% Compute energy spectrum
E_k = 0.5 * (u_k.*conj(u_k) + v_k.*conj(v_k));
E_K = zeros(n, 1);  % Correct size for energy bins

% Compute wave numbers
% k_x = zeros(size(1,n));
% k_y = zeros(size(1,n));
% dk = L/(2*n);

k_x =  fftfreq(n,(L/(n-1)));
k_y =  fftfreq(n,(L/(n-1)));

k_X = zeros(n,n);
k_Y = zeros(n,n);
for i = 1:n
k_X(i,:) = k_x;
k_Y(:,i) = k_y;
end

k_mag = sqrt(k_X.^2 + k_Y.^2);
% k_mag=zeros(n,n);

 for j = 1:n
    for i = 1:n
        k_id = round(k_mag(i,j))+1;
        if k_id <= n  % Prevent overflow
            E_K(k_id) = E_K(k_id) + E_k(i, j);  % Accumulate energy
        end
    end
 end






% Plot the energy spectrum
k = 1:n;  % Wavenumber range
plot(k, E_K, 'LineWidth', 1.5);
xlabel('Wavenumber k');
ylabel('Energy Spectrum E(k)');
title('Energy Spectrum');
grid on;

