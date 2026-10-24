%% Programmierprojekt 1 – Aufgabe 3: Kondition der Entzerrungsmatrix
% Zyklische Faltung y = C x mit zirkulanter Matrix C (erste Spalte = h, aufgefüllt).
% C = F^{-1} diag(h_hat) F  ->  Singulärwerte = |h_hat_k|  ->  kappa = max|h_hat| / min|h_hat|
% Die Funktion steht in kondition.m.

clear; close all; clc;

N = 16;
h = [1; 0.5; 0.25; -0.25];      % Echovektor aus Aufgabe 1b (siehe aufgabe1.m)


%% Kondition für h aus Aufgabe 1b + Kontrolle über die Matrix
kappa = kondition(h, N);

% Zirkulante Matrix explizit aufstellen: 1. Spalte = (h_0..h_{N-1}),
% 1. Zeile = (h_0, h_{N-1}, ..., h_1)
hp = [h; zeros(N - numel(h), 1)];
C  = toeplitz(hp, [hp(1); flipud(hp(2:end))]);

h_hat = fft(h, N);
F = fft(eye(N));                         % DFT-Matrix (F*x = fft(x))
fprintf('Kontrolle C = F^-1 diag(h_hat) F:   max. Abweichung %.1e\n', max(max(abs(F \ diag(h_hat) * F - C))));
fprintf('Kontrolle Singulärwerte = |h_hat|:  max. Abweichung %.1e\n', max(abs(sort(svd(C)) - sort(abs(h_hat)))));
fprintf('kappa über FFT (kondition.m) = %.4f\n', kappa);
fprintf('kappa über cond(C) (SVD)     = %.4f\n', cond(C));
[~, kmin] = min(abs(h_hat));  [~, kmax] = max(abs(h_hat));
fprintf('min |h_hat| = %.4f bei k = %d,  max |h_hat| = %.4f bei k = %d\n\n', ...
    abs(h_hat(kmin)), kmin - 1, abs(h_hat(kmax)), kmax - 1);


%% Familie h = (1, a), a in [-1, 0]
% h_hat_k = 1 + a * exp(-2*pi*i*k/N)
%  |h_hat_k| ist maximal bei k = N/2 (exp = -1):  1 - a = 1 + |a|
%  |h_hat_k| ist minimal bei k = 0   (exp = +1):  1 + a = 1 - |a|
%  -> kappa(a) = (1 - a) / (1 + a),  für a -> -1 geht kappa -> unendlich
a = linspace(-1, 0, 1001);
kap = arrayfun(@(ai) kondition([1; ai], N), a);
kap_formel = (1 - a) ./ (1 + a);
fprintf('Familie h = (1,a): max. rel. Abweichung Formel vs. numerisch: %.1e\n', ...
    max(abs(kap(2:end) - kap_formel(2:end)) ./ kap_formel(2:end)));
fprintf('kappa(a = -1) = %g  (C singulär: h_hat_0 = 1 - 1 = 0)\n', kap(1));
for ai = [-0.5 -0.9 -0.99 -0.999]
    fprintf('   a = %6.3f:  kappa = %8.1f\n', ai, kondition([1; ai], N));
end

figure('Name', 'Aufgabe 3 - Familie', 'Position', [100 100 1050 400]);
subplot(1, 2, 1);
semilogy(a, kap, 'LineWidth', 2); hold on;
semilogy(a(1:25:end), kap_formel(1:25:end), 'ro');
semilogy(-0.25, kondition(h, N), 'kp', 'MarkerSize', 12, 'MarkerFaceColor', 'y');
grid on; xlabel('a'); ylabel('\kappa_2(C)');
legend('numerisch (kondition.m)', '(1-a)/(1+a)', 'zum Vergleich: \kappa für h aus 1b', 'Location', 'northeast');
title('Kondition für h = (1, a)');
subplot(1, 2, 2);
kk = linspace(0, N, 400)';
farben = lines(4);
a_bsp = [0 -0.5 -0.9 -0.99];
for i = 1:numel(a_bsp)
    plot(kk, abs(1 + a_bsp(i) * exp(-2i * pi * kk / N)), '-', 'Color', farben(i, :)); hold on;
end
for i = 1:numel(a_bsp)
    plot(0:N-1, abs(fft([1; a_bsp(i)], N)), 'o', 'Color', farben(i, :), 'MarkerFaceColor', farben(i, :));
end
grid on; xlim([0 N]); xticks(0:15); xlabel('Unterträger k'); ylabel('|h\_hat_k|');
legend("a = " + string(a_bsp), 'Location', 'northwest');
title('Spektrum: Einbruch bei k = 0 für a \rightarrow -1');
speichere_plot('a3_familie.png');


%% Was bedeutet ein großes kappa bei Rauschen? (1) Datensatz 3
% Y_k = H_k X_k + E_k   (E_k: Rauschen auf Unterträger k)
% Entzerrt: Y_k / H_k = X_k + E_k / H_k  -> Fehler auf Unterträger k wird um 1/|H_k| verstärkt
H  = fft(h, N);
r3 = lese_datensatz('daten/datensatz3.csv');
y3 = r3(5:end);
Y3 = fft(y3);
X3 = Y3 ./ H;
X3_ent = qam_entscheiden(X3);
E = Y3 - H .* X3_ent;                   % Rauschen vor der Entzerrung
sigma = sqrt(mean(abs(y3 - ifft(H .* X3_ent)).^2));   % Rauschen im Zeitbereich (RMS)

fprintf('\nDatensatz 3: Rauschen vor Entzerrung |E_k| im Mittel %.3f\n', mean(abs(E)));
fprintf('             Fehler nach Entzerrung  |X_k - X_ent| = |E_k|/|H_k|\n');
R = corrcoef(abs(X3 - X3_ent), 1 ./ abs(H));
fprintf('             Korrelation Fehler <-> 1/|H_k|: %.2f\n', R(1, 2));
fprintf('             Rauschen im Zeitbereich: sigma = %.3f\n', sigma);

figure('Name', 'Aufgabe 3 - Datensatz 3', 'Position', [100 100 760 420]);
bar(0:N-1, [abs(E), abs(X3 - X3_ent)]); hold on;
plot(0:N-1, mean(abs(E)) ./ abs(H), 'k.-', 'MarkerSize', 14);
grid on; xticks(0:15); xlabel('Unterträger k'); ylabel('Betrag');
legend('|E_k| (Rauschen vor Entzerrung)', '|X''_k - Entscheidung| (nach Entzerrung)', ...
    'mittleres |E| / |H_k|', 'Location', 'northwest');
title(sprintf('Datensatz 3: Rauschverstärkung, wo |H_k| klein ist (\\kappa = %.2f)', kappa));
speichere_plot('a3_datensatz3_fehler.png');


%% (2) Simulation: gleiches Rauschen, verschiedene a
% Zufällige 16-QAM-Symbole über h = (1, a) senden, Rauschen wie in Datensatz 3,
% Zero-Forcing-Entzerrung, Fehler messen.
% Schranke aus der Vorlesung:  ||dx|| / ||x||  <=  kappa * ||dy|| / ||y||
rng(1);
M = 5000;                                   % Anzahl OFDM-Symbole
stufen = [-3 -1 1 3];
X = stufen(randi(4, N, M)) + 1i * stufen(randi(4, N, M));

a_test = [0 -0.25 -0.5 -0.75 -0.9 -0.95 -0.99];
fprintf('\nSimulation, %d OFDM-Symbole, sigma = %.3f:\n', M, sigma);
fprintf('%7s %8s %12s %12s %14s %12s %10s\n', 'a', 'kappa', 'rel.F. y', 'rel.F. x', ...
    'Verstärkung', 'max. Verst.', 'SER');
Xh_bsp = cell(size(a_test));
for i = 1:numel(a_test)
    Hi = fft([1; a_test(i)], N);
    y  = ifft(Hi .* X);                                     % Empfang ohne Rauschen
    dy = sigma / sqrt(2) * (randn(N, M) + 1i * randn(N, M));  % Rauschen
    Xh = fft(y + dy) ./ Hi;                                 % entzerren
    dx = ifft(Xh - X);                                      % Fehler im Zeitbereich

    rel_y = vecnorm(dy) ./ vecnorm(y);                      % pro OFDM-Symbol
    rel_x = vecnorm(dx) ./ vecnorm(ifft(X));
    ser   = mean(qam_entscheiden(Xh) ~= X, 'all');          % Symbolfehlerrate
    fprintf('%7.2f %8.1f %12.4f %12.4f %14.2f %12.2f %10.4f\n', a_test(i), kondition([1; a_test(i)], N), ...
        mean(rel_y), mean(rel_x), mean(rel_x ./ rel_y), max(rel_x ./ rel_y), ser);
    Xh_bsp{i} = Xh;
end
% "max. Verst." bleibt immer <= kappa (Schranke gilt), im Mittel ist die
% Verstärkung kleiner, weil nicht jedes Rauschen genau die schlechteste Richtung trifft.

figure('Name', 'Aufgabe 3 - Simulation', 'Position', [100 100 1100 380]);
zeige = [3 4 5];        % a = -0.5, -0.75, -0.9  (bei -0.99 fliegen die roten Punkte ganz aus dem Bild)
for s = 1:3
    subplot(1, 3, s);
    Xh = Xh_bsp{zeige(s)}(:, 1:300);
    plot(real(Xh(2:end, :)), imag(Xh(2:end, :)), '.', 'Color', [0.3 0.5 0.9], 'MarkerSize', 4); hold on;
    plot(real(Xh(1, :)), imag(Xh(1, :)), '.', 'Color', [0.9 0.2 0.2], 'MarkerSize', 6);
    [gi, gq] = meshgrid(stufen);
    plot(gi(:), gq(:), 'k+', 'MarkerSize', 10, 'LineWidth', 1.5);
    axis equal; axis([-6 6 -6 6]); grid on;
    title(sprintf('a = %.2f,  \\kappa = %.0f', a_test(zeige(s)), kondition([1; a_test(zeige(s))], N)));
    xlabel('Re'); ylabel('Im');
end
sgtitle('Entzerrte Symbole bei gleichem Rauschen (rot: Unterträger k = 0 mit kleinstem |h\_hat_k|)');
speichere_plot('a3_simulation.png');
