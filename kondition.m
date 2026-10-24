function kappa = kondition(h, N)
% KONDITION  Konditionszahl kappa_2(C) der zirkulanten Faltungsmatrix C
% zum Echovektor h (mit Nullen auf Länge N aufgefüllt).
%
% C = F^{-1} diag(h_hat) F  mit  h_hat = fft(h, N).  Da F/sqrt(N) unitär ist,
% sind die Singulärwerte von C genau |h_hat_k|, also
%   kappa_2(C) = max_k |h_hat_k| / min_k |h_hat_k|.
% Ist ein h_hat_k = 0, ist C singulär und kappa = Inf.

if nargin < 2
    N = 16;
end
h_hat = fft(h(:), N);
kappa = max(abs(h_hat)) / min(abs(h_hat));
end
