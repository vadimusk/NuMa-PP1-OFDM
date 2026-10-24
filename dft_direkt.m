function X = dft_direkt(x)
% DFT_DIREKT  Direkte Berechnung der DFT nach Definition:
%   X_k = sum_{n=0}^{N-1} x_n * exp(-2*pi*i*k*n/N),   k = 0..N-1
% Für jedes der N Ergebnisse wird eine Summe über N Terme gebildet -> O(N^2).
% Gleiche Konvention wie MATLABs fft.

x = x(:);
N = length(x);
n = (0:N-1)';
X = zeros(N, 1);
for k = 0:N-1
    X(k+1) = sum(x .* exp(-2i * pi * k * n / N));
end
end
