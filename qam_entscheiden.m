function Xd = qam_entscheiden(X)
% QAM_ENTSCHEIDEN  Rundet jeden Wert auf den nächsten 16-QAM-Punkt,
% d.h. Real- und Imaginärteil jeweils auf die nächste Zahl aus {-3,-1,1,3}.
% (Die Pilotwerte -3-3i, 1+1i, -1-1i, 3+3i liegen auch auf diesem Raster.)

stufe = @(v) min(max(2 * floor(v / 2) + 1, -3), 3);   % nächste ungerade Zahl, auf [-3,3] begrenzt
Xd = stufe(real(X)) + 1i * stufe(imag(X));
end
