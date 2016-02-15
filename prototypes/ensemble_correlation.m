function SPECTRAL_CORRELATION = ensemble_correlation(REGION_LIST, STEP)

% This function computes the sum of correlations for a list
% of 2-D arrays, and returns the result in the spectral Fourier domain.
%
% STEP is the amount by which to increment the list counter
% on each iteration. To correlate every region (i.e. time resolved
% data), set STEP = 1. For pair-wise (e.g., double-pulse) 
% correlations, set STEP = 2.

% Measure regions
[region_height, region_width, num_regions] = size(REGION_LIST);

% Allocate the correlation
SPECTRAL_CORRELATION = zeros(region_height, region_width);

% Compute the Fourier Transforms of all of the regions.
F = fft2(REGION_LIST);

% Region index
idx = 1;

% Do all the correlations
while idx < num_regions
        
    % Do the correlations between adjacent pairs, 
    % and add the result to the running sum of the
    % spectral correlation
    SPECTRAL_CORRELATION =  SPECTRAL_CORRELATION + ...
        F(:, :, idx) .* conj(F(:, :, idx + 1));

    % Increment the list counter
    idx = idx + STEP;
    
end


end