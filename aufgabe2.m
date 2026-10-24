%% Programmierprojekt 1 – Aufgabe 2: DFT vs. FFT – Laufzeitexperiment
% Wir messen die Laufzeit von
%   - dft_direkt   : eigene direkte DFT nach Definition, O(N^2)
%   - fft_rekursiv : eigene Radix-2-FFT (Cooley-Tukey), O(N log N)
%   - fft          : MATLAB-FFT (FFTW, kompilierter C-Code), O(N log N)
% für N = 2^4 ... 2^12 und bestimmen die Steigungen im log-log-Plot.
%
% Messmethode: siehe laufzeit.m (viele Wiederholungen, Minimum aus 5 Messungen).
% timeit haben wir zuerst probiert, das meldet bei fft aber "running too fast"
% und lieferte für kleine N unbrauchbare Werte (teilweise 0).

clear; close all; clc;

p  = 4:12;
Nw = 2.^p;
t_dft = zeros(size(Nw));
t_rek = zeros(size(Nw));
t_fft = zeros(size(Nw));

rng(0);
fprintf('%6s %14s %14s %14s %12s\n', 'N', 'DFT [s]', 'FFT eigen [s]', 'FFT MATLAB [s]', 'Fehler');
for j = 1:numel(Nw)
    N = Nw(j);
    x = randn(N, 1) + 1i * randn(N, 1);    % zufälliges komplexes Testsignal

    % Korrektheit prüfen (relativer Fehler gegenüber MATLAB-fft)
    fehler = max(norm(dft_direkt(x) - fft(x)), norm(fft_rekursiv(x) - fft(x))) / norm(fft(x));

    t_dft(j) = laufzeit(@() dft_direkt(x));
    t_rek(j) = laufzeit(@() fft_rekursiv(x));
    t_fft(j) = laufzeit(@() fft(x));
    fprintf('%6d %14.3e %14.3e %14.3e %12.1e\n', N, t_dft(j), t_rek(j), t_fft(j), fehler);
end


%% Empirische Steigungen
% Aus t ~ c * N^s folgt log t = log c + s * log N -> Gerade im log-log-Plot,
% Steigung s per linearer Regression (polyfit vom Grad 1).
% Bei kleinem N dominiert der feste Aufwand pro Aufruf (Funktionsaufruf,
% Speicher anlegen, ...), deshalb zusätzlich nur über N >= 2^8 gefittet.
gross = Nw >= 2^8;
steig = @(t, idx) polyfit(log10(Nw(idx)), log10(t(idx)), 1);

s = [steig(t_dft, true(size(Nw))); steig(t_rek, true(size(Nw))); steig(t_fft, true(size(Nw)))];
s_gross = [steig(t_dft, gross); steig(t_rek, gross); steig(t_fft, gross)];

fprintf('\nEmpirische Steigungen:\n');
fprintf('%-16s %12s %14s\n', '', 'alle N', 'nur N >= 256');
namen = {'DFT direkt', 'FFT eigen', 'FFT MATLAB'};
for i = 1:3
    fprintf('%-16s %12.2f %14.2f\n', namen{i}, s(i, 1), s_gross(i, 1));
end
% Theorie: DFT -> 2.  FFT: d log(N log N) / d log N = 1 + 1/ln(N)
fprintf('Theorie: DFT 2.00, FFT 1 + 1/ln(N) = %.2f (N=256) ... %.2f (N=4096)\n', ...
    1 + 1/log(256), 1 + 1/log(4096));


%% log-log-Plot
figure('Name', 'Aufgabe 2', 'Position', [100 100 760 520]);
loglog(Nw, t_dft, 'o-', 'LineWidth', 1.5); hold on;
loglog(Nw, t_rek, 's-', 'LineWidth', 1.5);
loglog(Nw, t_fft, 'd-', 'LineWidth', 1.5);
% Referenzgeraden, durch den letzten Messpunkt gelegt
loglog(Nw, t_dft(end) * (Nw / Nw(end)).^2, 'k--');
loglog(Nw, t_fft(end) * (Nw .* log2(Nw)) / (Nw(end) * log2(Nw(end))), 'k:');
grid on; xticks(Nw); xticklabels("2^{" + p + "}");
xlabel('N'); ylabel('Laufzeit pro Transformation [s]');
legend(sprintf('DFT direkt (Steigung %.2f)', s_gross(1, 1)), ...
       sprintf('FFT eigen, rekursiv (Steigung %.2f)', s_gross(2, 1)), ...
       sprintf('FFT MATLAB (Steigung %.2f)', s_gross(3, 1)), ...
       '\sim N^2', '\sim N log N', 'Location', 'northwest');
title('Laufzeit DFT vs. FFT (Steigungen gefittet für N \geq 2^8)');
speichere_plot('a2_laufzeit.png');


%% Abschätzung: OFDM-Symbole pro Sekunde
% Pro empfangenem OFDM-Symbol ist genau eine Transformation der Länge N nötig.
% -> maximal 1/t(N) Symbole pro Sekunde (nur FFT, ohne restliche Verarbeitung).
sys      = {'Projekt', 'LTE', '5G NR'};
N_sys    = [16 2048 4096];
benoetigt = [NaN 14000 28000];
N_genutzt = [16 1200 3276];

fprintf('\nOFDM-Symbole pro Sekunde (ein Prozessorkern, nur die Transformation):\n');
fprintf('%-8s %6s %10s %12s %12s %12s %16s\n', 'System', 'N', 'benötigt', ...
    'DFT', 'FFT eigen', 'FFT MATLAB', 'DFT nur genutzt');
for i = 1:3
    j = find(Nw == N_sys(i));
    % DFT nur für die genutzten Unterträger: N_genutzt Summen à N Terme
    dft_genutzt = 1 / (t_dft(j) * N_genutzt(i) / N_sys(i));
    fprintf('%-8s %6d %10g %12.0f %12.0f %12.0f %16.0f\n', sys{i}, N_sys(i), benoetigt(i), ...
        1 / t_dft(j), 1 / t_rek(j), 1 / t_fft(j), dft_genutzt);
end

fprintf('\nFaktor (geschafft / benötigt):\n');
for i = 2:3
    j = find(Nw == N_sys(i));
    fprintf('%-6s  DFT: %8.4f   FFT eigen: %8.2f   FFT MATLAB: %8.1f\n', sys{i}, ...
        (1 / t_dft(j)) / benoetigt(i), (1 / t_rek(j)) / benoetigt(i), (1 / t_fft(j)) / benoetigt(i));
end
