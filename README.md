# Programmierprojekt 1 – Entzerrung eines Funkkanals mit FFT

Numerische Mathematik für Elektrotechnik, WS 2026/27. MATLAB (getestet mit R2026a, keine Toolboxen nötig).

## Ausführen

MATLAB im Projektordner starten und nacheinander ausführen:

```matlab
aufgabe1   % Kanalschätzung und Entzerrung (Datensätze 1–3)
aufgabe2   % Laufzeitexperiment DFT vs. FFT (dauert ca. 20 s)
aufgabe3   % Kondition der Entzerrungsmatrix
```

Die Plots landen in `plots/`.

## Dateien

| Datei | Inhalt |
|---|---|
| `aufgabe1.m` | a) Pilotwerte prüfen, b) Echovektor h, c) Datensatz 3 entzerren, Zusatz: Präfix/ISI |
| `aufgabe2.m` | Laufzeitmessung, log-log-Plot, Steigungen, OFDM-Symbole/s vs. LTE/5G |
| `aufgabe3.m` | κ für h aus 1b, Familie h = (1, a), Auswirkung von Rauschen |
| `lese_datensatz.m` | CSV (Real-/Imaginärteil) → komplexer Spaltenvektor |
| `qam_entscheiden.m` | Rundet auf das nächste 16-QAM-Symbol |
| `dft_direkt.m` | DFT nach Definition, O(N²) |
| `fft_rekursiv.m` | Radix-2-FFT (Cooley-Tukey), O(N log N) |
| `kondition.m` | κ₂(C) = max\|ĥ_k\| / min\|ĥ_k\| |
| `laufzeit.m` | Zeitmessung mit Wiederholungen (Minimum aus 5 Messungen) |
| `speichere_plot.m` | Plot als PNG speichern |

## Ergebnisse

- **1a:** Nach Entfernen des Präfix und FFT sind genau k = 0, 4, 8, 12 belegt, mit −3−3i, 1+1i, −1−1i, 3+3i.
- **1b:** h = (1, 0.5, 0.25, −0.25). Die vier Pilot-Kanalwerte H₀, H₄, H₈, H₁₂ sind genau die DFT der Länge 4 von h → `h = ifft(H_pilot)`.
- **1c:** Entzerrt mit Xₖ = Yₖ / Hₖ, alle 16 Symbole liegen klar auf dem 16-QAM-Raster:
  `1+3i, 3−3i, −3+3i, −3−1i, 3−3i, 3−1i, 3+3i, −3+3i, −3+1i, −1−3i, 3+3i, 3+1i, −3−3i, 1−3i, −3−3i, −1+3i` (k = 0…15)
- **2:** Steigungen (N ≥ 256): DFT ≈ 1.9, MATLAB-FFT ≈ 1.13 (Theorie 1 + 1/ln N ≈ 1.12–1.18), eigene FFT ≈ 1.0 (Aufruf-Overhead).
  Mit direkter DFT schafft ein Kern nur ~16 Symbole/s bei N = 2048 (LTE braucht 14 000), mit `fft` ~170 000.
- **3:** κ(h aus 1b) = 3.34. Für h = (1, a) gilt κ = (1 − a)/(1 + a) → ∞ für a → −1 (bei a = −1 ist ĥ₀ = 0, C singulär).
  Rauschen wird auf Unterträgern mit kleinem |ĥₖ| um 1/|ĥₖ| verstärkt; relativ gilt ‖Δx‖/‖x‖ ≤ κ ‖Δy‖/‖y‖.
