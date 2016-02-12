function [x_grid_rect, y_grid_rect, ...
    x_grid_01, y_grid_01, ...
    x_grid_02, y_grid_02] = ...
    discreteWindowOffset(X_OLD, Y_OLD, U, V, ProcessingParameters) 
% Inputs: 
%   XGRID and YGRID are the matrices or vectors of the column and row 
%   grid points prior to shifting. These Should be in the format of meshgrid.
% 
%   U and V are the velocity fields by which to shift the grid points.
%   The size of U and V need not be the same as the size of XGRID and
%   YGRID. If the sizes of U and V are equal to the sizes of XGRID and
%   YGRID, then the grid will be shifted by exactly U and V; if those
%   sizes are not equal, then the shifts applied to XGRID and YGRID 
%   will be interpolated (linear) from U and V.  
%
% OUTPUTS
%   X_GRID_01 and Y_GRID_01 are the shifted grids for the first image in
%   the pair.
%   
%   X_GRID_02 and Y_GRID_02 are the shifted grids for the second image in
%   the pair.
%   

% Read image dimensions
imageHeight = ProcessingParameters.Images.Height;
imageWidth  = ProcessingParameters.Images.Width;

% Grid spacing parameters
gridSpacingY = ProcessingParameters.Grid.Spacing.Y;
gridSpacingX = ProcessingParameters.Grid.Spacing.X;

% Make sure the grid buffer is at least half the size of the interrogation region
gridBufferY = ProcessingParameters.Grid.Buffer.Y;
gridBufferX = ProcessingParameters.Grid.Buffer.X;

% Create the initial rectangular grid
[x_grid_rect, y_grid_rect] = ...
    gridImage([imageHeight, imageWidth], ...
    [gridSpacingY gridSpacingX], gridBufferY, gridBufferX);

% Interpolate the shift-field onto the input grid. U and V are the
% non-integer horizontal and vertical displacements by whose rounded values
% the grid is shifted.
% The first  "linear" is the interpolation method and
% The second "linear" is the extrapolation method.
% We are using scatteredInterpolant because interp2 doesn't have a
% good extrapolation method input. 
interpolant_U = scatteredInterpolant(Y_OLD,X_OLD, U, 'linear','linear');
interpolant_V = scatteredInterpolant(Y_OLD,X_OLD, V, 'linear','linear');
gridShiftX = interpolant_U(y_grid_rect, x_grid_rect);
gridShiftY = interpolant_V(y_grid_rect, x_grid_rect);

% Shift the grid points from the first image by -1/2 times the input
% displacement field.
x_grid_01 = x_grid_rect - round(1/2 * gridShiftX);
y_grid_01 = y_grid_rect - round(1/2 * gridShiftY);

% Shift the grid points from the second image by +1/2 times the input
% displacement field
x_grid_02 = x_grid_rect + round(1/2 * gridShiftX);
y_grid_02 = y_grid_rect + round(1/2 * gridShiftY);


end


