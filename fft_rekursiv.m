function alpha = fft_rekursiv(f)
% FFT_REKURSIV  Schnelle Fourier-Transformation nach Algorithmus 1.1 im Skript
% (Lemma 1.21), rekursiv, O(N log N). Berechnet alpha = F_N * f (wie fft).
% Die Länge N von f muss eine Zweierpotenz sein (N = 2^p, p >= 1).
%
% Idee (Lemma 1.21): eine DFT der Länge N -> zwei DFTs der Länge N/2
%   gerade Einträge:   alpha_{2m}   = DFT_{N/2}( f_l + f_{l+N/2} )
%   ungerade Einträge: alpha_{2m+1} = DFT_{N/2}( (f_l - f_{l+N/2}) * w^l ),  w = exp(-2*pi*i/N)

f = f(:);
N = length(f);
if N == 2
    alpha = [f(1) + f(2); f(1) - f(2)];     % F_2 = [1 1; 1 -1]
    return
end

w = exp(-2i * pi * (0:N/2-1)' / N);         % Potenzen der Einheitswurzel omega_N
f_plus  = f(1:N/2) + f(N/2+1:N);
f_minus = (f(1:N/2) - f(N/2+1:N)) .* w;

alpha = zeros(N, 1);
alpha(1:2:N) = fft_rekursiv(f_plus);         % alpha_0, alpha_2, alpha_4, ...
alpha(2:2:N) = fft_rekursiv(f_minus);        % alpha_1, alpha_3, alpha_5, ...
end
