function X = fft_rekursiv(x)
% FFT_REKURSIV  Radix-2-FFT nach Cooley-Tukey (rekursiv), O(N log N).
% Die Länge von x muss eine Zweierpotenz sein.
%
% Idee: DFT der Länge N = zwei DFTs der Länge N/2 (gerade / ungerade Indizes)
%   X_k       = G_k + w^k * U_k
%   X_{k+N/2} = G_k - w^k * U_k,     w = exp(-2*pi*i/N),  k = 0..N/2-1

x = x(:);
N = length(x);
if N == 1
    X = x;          % DFT der Länge 1 ist der Wert selbst
    return
end

G = fft_rekursiv(x(1:2:end));      % x_0, x_2, x_4, ...
U = fft_rekursiv(x(2:2:end));      % x_1, x_3, x_5, ...
w = exp(-2i * pi * (0:N/2-1)' / N); % Twiddle-Faktoren
X = [G + w .* U; G - w .* U];
end
