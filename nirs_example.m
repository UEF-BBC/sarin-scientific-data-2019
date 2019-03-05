% A simple example script for the equine dataset.
% This script demonstrates how to utilize the NIRS and reference data to build 
% a simple partial least squares model.
%
% Mithilesh Prakash, Jari Torniainen
% Department of Applied Physics, Univeristy of Eastern Finland, 2019

clear all; close all; clc;

%% Load NIR and reference data
data_structure = load('nirs_and_references.mat');
dataset = data_structure.dataset;

wavelength = dataset(1).wavelength;

% Average the spectra from each measurement point (3 spectra per 869 locations)
for measurement = 1:length(dataset)
    spectra_raw = dataset(measurement).raw_spectra;
    spectra_raw_averaged(measurement, :) = mean(spectra_raw);
end

%load a reference parameter: e.g. Cartilage thickness
reference_value = [dataset.cartilage_thickness].';

% Remove instances with missing reference data
spectra_raw_averaged = spectra_raw_averaged(find(~isnan(reference_value)), :);
reference_value = reference_value(find(~isnan(reference_value)), :);

%% Preprocess the spectral data
Preprocess = 3;
window_length = 41;
polynomial_order = 3;

[~, g] = sgolay(polynomial_order, window_length);
half_window = (window_length + 1) / 2 - 1;
for ii = 1:size(spectra_raw_averaged, 1)
    for n = (window_length + 1) / 2:size(spectra_raw_averaged, 2) - (window_length + 1) / 2
        spectra_preprocess(ii, n - (window_length + 1) / 2 + 1) = dot(g(:, Preprocess), spectra_raw_averaged(ii, n - half_window:n + half_window));
    end
end
wavelength = wavelength((window_length + 1)./2:end - (window_length + 1)./2);

%% Simple Partial least squares regression model (with 1-10 components)
maximum_num_components = 10;
num_of_cv_folds = 10;
[~, ~, ~, ~, beta_values, pct_var, msep, stats] = plsregress(spectra_preprocess, reference_value, maximum_num_components, 'CV', num_of_cv_folds);

figure();

subplot(2, 2, 1);
plot(1:10, cumsum(100 * pct_var(2, :)), '-bo');
xlabel('Number of PLS components');
ylabel('Variance explained in reference (%)');
title('Variance Plot');

subplot(2, 2, 2) ;
plot(0:10, msep(2, :));
xlabel('Number of components');
ylabel('Mean Squared Prediction Error');
title('MSE Plot');
legend({'PLSR'}, 'location', 'NE');

subplot(2, 2, 3);
plot(wavelength, stats.W(:, 1:3), '-');
xlabel('Wavelengths');
ylabel('PLS Weights');
title('Weights on Wavelengths (three component PLS model)');
legend({'1st Component', '2nd Component', '3rd Component'}, 'location', 'NW');

%% Computing the fitted response and residuals.
reference_fitted = [ones(size(spectra_preprocess, 1), 1) spectra_preprocess] * beta_values;
residuals = reference_value - reference_fitted;
subplot(2, 2, 4);
stem(residuals);
xlabel('Observation');
ylabel('Residual');
title('Residual Plot');
