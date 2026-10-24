%% Programmierprojekt 1 – Aufgabe 1: Kanalschätzung und Entzerrung
% Numerische Mathematik für Elektrotechnik, WS 2026/27
%
%  a) Datensatz 1: Präfix entfernen, FFT, Pilotwerte prüfen
%  b) Datensatz 2: Echovektor h mit Pilotwerten + Faltungssatz bestimmen
%  c) Datensatz 3: mit h entzerren und die 16 Symbole plotten
%
% Konvention: Unterträger k = 0..15 liegt in MATLAB bei Index k+1.

clear; close all; clc;

N  = 16;            % Anzahl Unterträger (= FFT-Länge)
Lp = 4;             % Länge des zyklischen Präfix
k  = (0:N-1)';

k_pilot = [0 4 8 12];
ip      = k_pilot + 1;                               % MATLAB-Indizes der Piloten
X_pilot = lese_datensatz('daten/pilotsymbole.csv');  % 16 Werte, außer bei k_pilot alle 0


%% a) Pilotwerte prüfen (Datensatz 1)
r1 = lese_datensatz('daten/datensatz1.csv');    % 20 Abtastwerte (Zeitbereich)

% Der zyklische Präfix ist eine Kopie der letzten 4 Werte des Symbols
fprintf('a) Präfix = letzte %d Werte des Symbols?  max. Abweichung: %.1e\n', ...
    Lp, max(abs(r1(1:Lp) - r1(end-Lp+1:end))));

x1 = r1(Lp+1:end);  % Präfix entfernen -> 16 Werte
X1 = fft(x1);       % in den Frequenzbereich: X1(k+1) = Wert auf Unterträger k

disp(table(k, round(X1, 12), X_pilot, 'VariableNames', {'k', 'FFT_Datensatz1', 'Pilotwerte'}))
fprintf('   max |FFT - Pilotwerte| = %.1e\n', max(abs(X1 - X_pilot)));
fprintf('   belegte Unterträger (|X_k| > 1e-10): k = %s\n\n', mat2str(k(abs(X1) > 1e-10)'));

figure('Name', 'Aufgabe 1a');
stem(k, abs(X1), 'filled'); hold on;
plot(k_pilot, abs(X_pilot(ip)), 'ro', 'MarkerSize', 12, 'LineWidth', 1.5);
for j = ip
    text(k(j) + 0.3, abs(X1(j)) + 0.2, sprintf('%.0f%+.0fi', real(X1(j)), imag(X1(j))));
end
grid on; xlim([-0.5 15.5]); ylim([0 5]); xticks(0:15);
xlabel('Unterträger k'); ylabel('|X_k|');
legend('|FFT(Datensatz 1)|', 'bekannte Pilotwerte');
title('a) Datensatz 1 im Frequenzbereich');
speichere_plot('a1a_pilotwerte.png');


%% b) Echovektor bestimmen (Datensatz 2)
r2 = lese_datensatz('daten/datensatz2.csv');
y2 = r2(Lp+1:end);
Y2 = fft(y2);

% Faltungssatz:  y = h (*) x  (zyklische Faltung)   <=>   Y_k = H_k * X_k
% mit H = fft(h, N). Auf den Pilot-Unterträgern kennen wir X_k und Y_k:
H_pilot = Y2(ip) ./ X_pilot(ip);

% Auf den leeren Unterträgern ist Y_k = H_k * 0 = 0 -> dort erfährt man nichts über H_k
fprintf('b) max |Y_k| auf leeren Unterträgern: %.1e\n', max(abs(Y2(setdiff(1:N, ip)))));
disp(table(k_pilot', Y2(ip), X_pilot(ip), H_pilot, 'VariableNames', {'k', 'Y_k', 'X_k', 'H_k'}))

% h hat nur Länge 4:  H_k = sum_{n=0}^{3} h_n * exp(-2*pi*i*k*n/16).
% Für k = 4m (m = 0..3) ist exp(-2*pi*i*4m*n/16) = exp(-2*pi*i*m*n/4),
% also ist (H_0, H_4, H_8, H_12) genau die DFT der Länge 4 von h.
% -> h ist die inverse DFT der Länge 4 der vier Pilot-Kanalwerte:
h = ifft(H_pilot);

% Kontrolle: dasselbe als lineares 4x4-Gleichungssystem A*h = H_pilot
n = 0:3;
A = exp(-2i * pi * k_pilot' * n / N);
h_lgs = A \ H_pilot;
fprintf('   Abweichung ifft vs. LGS-Lösung: %.1e\n', max(abs(h - h_lgs)));
fprintf('   max |Im(h)| = %.1e  -> h ist reell\n', max(abs(imag(h))));
h = real(h);
fprintf('   h = [%s]\n', strjoin(compose('%g', h'), ', '));

H = fft(h, N);      % Kanal jetzt auf ALLEN 16 Unterträgern bekannt

% Probe ohne FFT, direkt über die Definition der zyklischen Faltung:
%   y_j = sum_n h_n * x_{(j-n) mod N}     (x = gesendetes Pilotsymbol = x1 aus a)
h16 = [h; zeros(N - 4, 1)];
y2_probe = zeros(N, 1);
for j = 0:N-1
    for m = 0:N-1
        y2_probe(j+1) = y2_probe(j+1) + h16(m+1) * x1(mod(j - m, N) + 1);
    end
end
fprintf('   Probe: max |h (*) x1 - y2| = %.1e\n\n', max(abs(y2_probe - y2)));

figure('Name', 'Aufgabe 1b', 'Position', [100 100 1000 380]);
subplot(1, 2, 1);
stem(0:3, h, 'filled', 'LineWidth', 1.5);
grid on; xlim([-0.5 3.5]); xticks(0:3); ylim([-0.5 1.2]);
xlabel('n (Verzögerung in Abtastwerten)'); ylabel('h_n'); title('Echovektor h');
subplot(1, 2, 2);
kk = linspace(0, N, 500)';
plot(kk, abs(exp(-2i * pi * kk * (0:3) / N) * h), '-', 'Color', [0.6 0.6 0.6]); hold on;
plot(k, abs(H), 'bo', 'MarkerFaceColor', 'b');
plot(k_pilot, abs(H_pilot), 'rs', 'MarkerSize', 12, 'LineWidth', 1.5);
grid on; xlim([0 N]); xticks(0:15);
xlabel('Unterträger k'); ylabel('|H_k|');
legend('|H| kontinuierlich', '|H_k| = |fft(h,16)|', 'aus Piloten gemessen', 'Location', 'south');
title('Kanal im Frequenzbereich');
speichere_plot('a1b_echovektor.png');


%% c) Unbekannte Werte bestimmen (Datensatz 3)
r3 = lese_datensatz('daten/datensatz3.csv');
y3 = r3(Lp+1:end);
Y3 = fft(y3);

X3 = Y3 ./ H;               % Entzerrung: X_k = Y_k / H_k   (entspricht x = C^{-1} y)
X3_ent = qam_entscheiden(X3);   % auf das 16-QAM-Raster runden

disp('c) Entzerrte Symbole aus Datensatz 3:')
disp(table(k, round(Y3, 3), round(X3, 3), X3_ent, round(abs(X3 - X3_ent), 3), ...
    'VariableNames', {'k', 'Y_k', 'X_k_entzerrt', 'Entscheidung', 'Abstand'}))

figure('Name', 'Aufgabe 1c', 'Position', [100 100 1000 460]);
[gi, gq] = meshgrid([-3 -1 1 3]);
subplot(1, 2, 1);
plot(gi(:), gq(:), 'k+', 'MarkerSize', 10); hold on;
plot(real(Y3), imag(Y3), 'ro', 'MarkerFaceColor', 'r');
text(real(Y3) + 0.12, imag(Y3) + 0.18, string(k), 'FontSize', 8);
axis equal; grid on; axis([-6.5 6.5 -6.5 6.5]);
xlabel('Re'); ylabel('Im'); title('vor der Entzerrung: Y_k');
subplot(1, 2, 2);
plot(gi(:), gq(:), 'k+', 'MarkerSize', 10); hold on;
plot(real(X3), imag(X3), 'bo', 'MarkerFaceColor', 'b');
text(real(X3) + 0.12, imag(X3) + 0.18, string(k), 'FontSize', 8);
axis equal; grid on; axis([-4.5 4.5 -4.5 4.5]);
xlabel('Re'); ylabel('Im'); title('nach der Entzerrung: X''_k = Y_k / H_k');
sgtitle('c) Die 16 übertragenen Symbole (Zahl = Unterträger k)');
speichere_plot('a1c_symbole.png');


%% Zusatz: Warum den Präfix wegwerfen? (Datensatz 2 und 3 direkt nacheinander)
% Physikalisch ist der Kanal eine LINEARE Faltung des gesamten Sendesignals.
% Das Echo von Symbol 2 "läuft" in den Anfang von Symbol 3 hinein
% (Intersymbolinterferenz). Wir bauen das nach und vergleichen mit Datensatz 3.
s2 = [x1(end-Lp+1:end); x1];       % Symbol 2 = Pilotsymbol, mit Präfix
x3 = ifft(X3_ent);                 % gesendetes Symbol 3 (aus den entschiedenen Werten)
s3 = [x3(end-Lp+1:end); x3];       % Symbol 3 mit Präfix

r_mit  = conv([s2; s3], h);  r_mit  = r_mit(21:40);   % mit Vorgänger-Symbol
r_ohne = conv(s3, h);        r_ohne = r_ohne(1:20);   % ohne Vorgänger-Symbol

fprintf('Zusatz: |Simulation - Datensatz 3| pro Abtastwert (Präfix = Werte 1..4)\n');
disp(table((1:20)', round(abs(r_mit - r3), 3), round(abs(r_ohne - r3), 3), ...
    'VariableNames', {'Abtastwert', 'mit_Symbol_2_davor', 'ohne_Symbol_2'}))
% -> ohne Vorgänger passen nur die ersten 3 Werte (alle im Präfix) nicht.
%    Nach dem Abschneiden des Präfix ist der Rest exakt eine ZYKLISCHE Faltung.

sigma = sqrt(mean(abs(y3 - ifft(H .* X3_ent)).^2));
fprintf('geschätzte Rauschstärke in Datensatz 3 (Zeitbereich, RMS): %.3f\n', sigma);
