% An  example script for the equine dataset.
% This script demonstrates how to compute the 2d maps of amide I and proteoglycan (PG) 
% from the raw FTIR data.
%
% Jari Torniainen, Lassi Rieppo, Jaakko Sarin
% Department of Applied Physics, Univeristy of Eastern Finland, 2019

clear all; close all; clc;

load('ftir_raw.mat');  % Load FTIR measurements. Using only the first sample of the dataset.
load('horse_reference_blc.mat');  % Load reference spectra.

% Unpack data and reshape it to row vectors.
wave_raw = ftir_raw(1).wave;
spectra_raw = ftir_raw(1).data;
spectra_raw_ = reshape(spectra_raw, size(spectra_raw, 1) * size(spectra_raw, 2), size(spectra_raw, 3));

% Preprocess data. Preprocessing is done using:
% "Resonant Mie Scattering EMSC (RMieS-EMSC) correction for infrared spectroscopy"
% https://github.com/GardnerLabUoM/RMieS
% Make sure to download the toolbox and add it to the path.
RMIES_PATH = 'RMieS';
addpath(RMIES_PATH);
correction_options = [ ...
                        0    ;      % 1. Desired resolution, (0 keeps original resolution)
                        700  ;      % 2. Lower wavenumber range (min value is 1000)
                        4000 ;      % 3. Upper wavenumber range (max value is 4000)
                        5    ;      % 4. Number of iterations
                        2    ;      % 5. Mie theory option (smooth or RMieS)
                        8    ;      % 6. Number of principal components used (8 recommended)
                        2    ;      % 7. Lower range for scattering particle diameter / um
                        8    ;      % 8. Upper range for scattering particle diamter / um
                        1.1  ;      % 9. Lower range for average refractive index
                        1.5  ;      % 10. Upper range for average refractive index
                        10   ;      % 11. Number of values for each scattering parameter default 10
                        1    ;      % 12. Orthogonalisation, 0 = no, 1 = yes. (1 recommended)
                        0   ];      % 13. Which reference spectrum, 1 = Matrigel, 2 = Simulated, 0 = custom

[wave, spectra, history] = RMieS_EMSC_v5(wave_raw, spectra_raw_, correction_options, WN_ref, horse_reference_blc);

% Find wavenumber regions for amide I and PG
AMIDE_I_BAND = [1590, 1720];
PG_BAND = [984, 1140];
amide_i_idx = find_band(wave, AMIDE_I_BAND);
pg_idx = find_band(wave, PG_BAND);

% Find amide I and pg values by evaluating the integral across the two regions.
amide_peaks = trapz(wave(amide_i_idx), spectra(:, amide_i_idx)');
pg_peaks = trapz(wave(pg_idx), spectra(:, pg_idx)');

% Remove saturated and misaligned spectra
CUTOFF_HIGH = 2.5; % Threshold for saturated amide I peak
CUTOFF_LOW = 0.05; % Threshold for non-existant amide I peak
bad_data = sum(spectra(:, amide_i_idx) > CUTOFF_HIGH, 2) > 0 | max(spectra(:, amide_i_idx), [], 2) < CUTOFF_LOW;
amide_peaks(bad_data) = NaN
pg_peaks(bad_data) = NaN

% Reshape amide I and pg peak values to 2d maps.
amide_map = reshape(amide_peaks, size(spectra_raw, 1), size(spectra_raw, 2));
pg_map = reshape(pg_peaks, size(spectra_raw, 1), size(spectra_raw, 2));

% Visualize results.
figure();
subplot(2, 2, 1);  % Plot raw spectra and highlight rejected spectra (in red).
hold on;
plot(wave_raw, spectra_raw_', '-', 'color', [.6, .6, .6, .05]);
plot(wave_raw, mean(spectra_raw_, 1), '-k', 'linewidth', 2);
plot(wave_raw, spectra_raw_(bad_data, :)', '-r', 'linewidth', 1);
hold off;
axis('tight');
title('Raw FTIR spectra');
xlabel('Wavenumber [cm^{-1}]');
ylabel('Absorption [A.U.]');
grid();

subplot(2, 2, 2);  % Plot corrected spectra and highlight regions for amide I and pg.
hold on;
plot(wave, spectra', '-', 'color', [.6, .6, .6, .05]);
plot(wave, mean(spectra, 1), '-k', 'linewidth', 2);
axis('tight');
ylims = get(gca, 'ylim');
patch([AMIDE_I_BAND, fliplr(AMIDE_I_BAND)], [ylims(1), ylims(1), ylims(2), ylims(2)], [1, .4, .4], 'facealpha', .3, 'edgecolor', [1, .4, .4]);
patch([PG_BAND, fliplr(PG_BAND)], [ylims(1), ylims(1), ylims(2), ylims(2)], [.4, .4, 1], 'facealpha', .3, 'edgecolor', [.4, .4, 1]);
text(mean(AMIDE_I_BAND), ylims(2), 'Amide I', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontweight', 'bold', 'color', 'r');
text(mean(PG_BAND), ylims(2), 'PG', 'horizontalalignment', 'center', 'verticalalignment', 'top', 'fontweight', 'bold', 'color', 'b');
hold off;
title('Corrected FTIR spectra');
xlabel('Wavenumber [cm^{-1}]');
ylabel('Absorption [A.U.]');
grid();

subplot(2, 2, 3);  % Plot the resulting 2D amide I map.
imagesc(amide_map');
colorbar;
title('Amide I map');
xlabel('Horizontal position');
ylabel('Vertical position');

subplot(2, 2, 4);  % Plot the resulting 2D PG map.
imagesc(pg_map');
colorbar;
title('PG map');
xlabel('Horizontal position');
ylabel('Vertical position');

function mask = find_band(wave, band_limits)
    % Returns a (binary) mask for selecting the wavenumber band specified in band_limits
    band_limits = sort(band_limits);
    mask = wave >= band_limits(1) & wave <= band_limits(2);
end
