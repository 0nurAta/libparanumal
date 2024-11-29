% Flow variables
L  = 2*pi;    % Domain Length
nu = 0.000625;  % Viscosity
t  = 10;       % Time
n = 24;
% Calculate the exact velocity
x = linspace(0,L,n);  % Corrected range
y = linspace(0,L,n);
n = size(x,2);
uexact = zeros(n, n);
vexact = zeros(n, n);

for i = 1:n
    for j = 1:n
        uexact(i, j) = sin(x(i)) * cos(y(j)) * exp(-2 * nu * t);
        vexact(i, j) = -cos(x(i)) * sin(y(j)) * exp(-2 * nu * t);
    end
end

% Apply Fourier Transform
% u_k = abs(fft2(uexact)) / (n * n);  % Normalize
% v_k = abs(fft2(vexact)) / (n * n);  % Normalize

u_k = fft(uexact,n,1);  % Normalize
u_k = fft(u_k,n,2) / (n * n);
v_k = fft(vexact,n,1);  % Normalize
v_k = fft(v_k,n,2) / (n * n);

% u_k = fft2(uexact)/ (n * n);  % Normalize
% v_k = fft2(vexact)/ (n * n);  % Normalize


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

