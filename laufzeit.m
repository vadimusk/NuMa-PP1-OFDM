function t = laufzeit(f)
% LAUFZEIT  Laufzeit eines einzelnen Aufrufs von f (Funktionshandle) in Sekunden.
%
% Sehr schnelle Aufrufe (z.B. fft mit N = 16, < 1 Mikrosekunde) kann tic/toc
% nicht einzeln messen. Deshalb: f so oft hintereinander ausführen, bis
% mindestens 50 ms zusammenkommen, und durch die Anzahl teilen.
% Das Ganze 5-mal; genommen wird das Minimum (am wenigsten durch andere
% Prozesse gestört).

f();                        % Aufwärmen (JIT-Compiler, Speicher, FFTW-Plan)

R = 1;                      % Anzahl Wiederholungen pro Messung bestimmen
while true
    tic;
    for r = 1:R, f(); end
    if toc >= 0.05, break; end
    R = 2 * R;
end

t = Inf;
for v = 1:5
    tic;
    for r = 1:R, f(); end
    t = min(t, toc / R);
end
end
