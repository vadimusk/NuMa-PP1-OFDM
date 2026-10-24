function z = lese_datensatz(dateiname)
% LESE_DATENSATZ  Liest eine Zwei-Spalten-Textdatei (Realteil Imaginärteil)
% und gibt die Werte als komplexen Spaltenvektor zurück.
% Die erste Zeile ist ein Kommentar (#) und wird übersprungen.

M = readmatrix(dateiname, 'FileType', 'text', 'CommentStyle', '#', 'Delimiter', ' ');
z = M(:, 1) + 1i * M(:, 2);
end
