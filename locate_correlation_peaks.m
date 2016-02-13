%#codegen
function PEAK_MATRIX = locate_correlation_peaks(SPATIAL_CORRELATION)
% This function accepts a spatial correlation plane
% and returns the plane with peaks extracted. 
% This tiny function is separate from the main subpixel code
% so that it can be compiled easily.
% 
% Compile with: codegen -config coder.config('mex') locate_correlation_peaks -args {coder.typeof(1.00, [inf, inf])};
%
% INPUTS:
%	SPATIAL_CORRELATION = [M x N] matrix containing a spatial correlation plane
%
% OUTPUTS:
%	PEAK_MATRIX = [M x N] matrix containing the locations and magnitudes of the identified peaks.

% Check if the correlation was identically zero
if abs(max(SPATIAL_CORRELATION(:))) ~= 0
    % Locate peaks using imregionalmax
    PEAK_MATRIX = imregionalmax(SPATIAL_CORRELATION) .* SPATIAL_CORRELATION;
else
    PEAK_MATRIX = zeros(size(SPATIAL_CORRELATION));   
end

end
