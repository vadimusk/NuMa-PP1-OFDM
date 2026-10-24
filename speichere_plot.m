function speichere_plot(dateiname)
% SPEICHERE_PLOT  Aktuelle Abbildung als PNG nach plots/ speichern.
% Immer im hellen Design (sonst übernimmt MATLAB ab R2025a den Dark Mode des Systems).

if ~exist('plots', 'dir'), mkdir('plots'); end
if ~isMATLABReleaseOlderThan("R2025a")
    theme(gcf, 'light');
end
exportgraphics(gcf, fullfile('plots', dateiname), 'Resolution', 150);
end
