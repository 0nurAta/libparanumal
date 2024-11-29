% Solution Information
L  = 2*pi;

% Import Solution Data
u_soln = readtable("deneme.csv",Range="A2:A185194",ReadVariableNames=false);
v_soln = readtable("deneme.csv",Range="B2:B185194",ReadVariableNames=false);
w_soln = readtable("deneme.csv",Range="C2:C185194",ReadVariableNames=false);
u_soln = table2array(u_soln(:,1));
v_soln = table2array(v_soln(:,1));
w_soln = table2array(w_soln(:,1));

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

u_k = fft(u,n,1);  % Normalize
u_k = fft(u_k,n,2);
u_k = fft(u_k,n,3) / (n*n*n);
v_k = fft(v,n,1);  % Normalize
v_k = fft(v_k,n,2);
v_k = fft(v_k,n,3) / (n*n*n);
w_k = fft(w,n,1);  % Normalize
w_k = fft(w_k,n,2);
w_k = fft(w_k,n,3) / (n*n*n);


% u_k = fftshift(u_k);
% v_k = fftshift(v_k);
% Compute energy spectrum
E_k = 0.5 * (u_k.*conj(u_k) + v_k.*conj(v_k)+ w_k.*conj(w_k));
E_K = zeros(n, 1);  % Correct size for energy bins

% Compute wave numbers
% k_x = zeros(size(1,n));
% k_y = zeros(size(1,n));
% dk = L/(2*n);

k_x =  fftfreq(n,(L/(n-1)));
k_y =  fftfreq(n,(L/(n-1)));
k_z =  fftfreq(n,(L/(n-1)));

k_X = zeros(n,n);
k_Y = zeros(n,n);
k_Z = zeros(n,n);
for i = 1:n
k_X(i,:) = k_x;
k_Y(i,:) = k_y;
k_Z(i,:) = k_z;
end

k_X3 = zeros(n,n,n);
k_Y3 = zeros(n,n,n);
k_Z3 = zeros(n,n,n);

for i = 1:n
k_X3(i,:,:) = k_X;
k_Y3(i,:,:) = k_Y;
k_Z3(i,:,:) = k_Z;
end

k_mag = sqrt(k_X3.^2 + k_Y3.^2+ k_Z3.^2);
% k_mag=zeros(n,n);
kMag = zeros(n);
 for j = 1:n
    for i = 1:n
        for l = 1:n
        k_id = round(k_mag(i,j,l))+1;
            if k_id <= n  % Prevent overflow
            kMag(k_id) = round(k_mag(i,j,l))+1;
            E_K(k_id) = E_K(k_id) + E_k(i, j,l);  % Accumulate enery
            end
        end
    end
 end

% Plot the energy spectrum
k = 1:n;  % Wavenumber range
plot(k, E_K, 'LineWidth', 1.5);
loglog(E_K)
xlabel('Wavenumber k');
ylabel('Energy Spectrum E(k)');
title('Energy Spectrum');
grid on;

