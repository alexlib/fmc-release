function pivFullField(FILEPATHS, JOBFILE)

% Check for whether to use compiled codes
% Only use compiled codes on Linux.
% This can be changed if you want to re-compile the codes
% on your own platform, i.e. if you figure out
% how to compile mex files on Mac, etc.
COMPILED = JOBFILE.JobOptions.RunCompiled;

% Number of PIV passes
% Take the minimum of the requested number of passes and the number of pass
% parameter structures specified.
number_of_passes = min(JOBFILE.JobOptions.NumberOfPasses, length(JOBFILE.Parameters.Processing));

% Truncate the jobfile processing parameters
% down to the number of passes specified by
% the variable "numberOfPasses". 
JOBFILE.Parameters.Processing = JOBFILE.Parameters.Processing(1 : number_of_passes);

% Read in the clean images.
image1_import = double(imread(FILEPATHS.FirstImagePath));
image2_import = double(imread(FILEPATHS.SecondImagePath));

% Image dimensions
[imageHeight, imageWidth, num_channels] = size(image1_import);

% Extract the channel
if num_channels > 1   
    % Make sure the color channel is specified
    if isfield(JOBFILE.Parameters.Images, 'ColorChannel')
        % Read the requested color channel
        channel = JOBFILE.Parameters.Images.ColorChannel;
        % Extract the color channel from the image.
        image1_raw = image1_import(:, :, channel);
        image2_raw = image2_import(:, :, channel);
    else
        % Default to taking the first channel if none is specified.
        image1_raw = image1_import(:, :, 1);
        image2_raw = image2_import(:, :, 1);
    end
    
else  % else (if num_channels > 1)
    
    % If there's only one channel
    % then set "image1_raw" to be the whole image, etc.
    image1_raw = image1_import;
    image2_raw = image2_import;
    
end

% Check whether to re-start from a previously existing velocity field.
startFromExistingField = JOBFILE.JobOptions.StartFromExistingField;

% If starting from a previously existing velocity field,
% load that velocity field.
if startFromExistingField && exist(FILEPATHS.OutputFilePath, 'file') 
    
    % Load the existing velocity field if it exists.
    load(FILEPATHS.OutputFilePath);
        
    % Set the current pass number to the first specified pass
    % If "start from existing" path is specified and 
    % a starting-pass number is specified less than 2,
    % then default to starting at pass 2. This 
    % deals with inputs of 0, negative numbers, 
    % or confusion about what "starting pass" means.
    % Presumably if the user specifies starting from an existing
    % pass, they want to pick things up starting with pass 2
    % (starting on pass 1 would just repeat pass 1, which is the same
    % thing as startFromExistingField = 0)
    if JOBFILE.JobOptions.StartPass < 2
        piv_pass_number = 2;
    else
        piv_pass_number = JOBFILE.JobOptions.StartPass;
    end
    
    % This sets a loop counter for the PIV passes.
    p = piv_pass_number - 1;
    
    % Create the function variables from the saved data
    % This just entails flipping all the coordinate and 
    % displacement data
    % I'm not sure why all this flipping is done...
    for pass = 1 : length(X)
        gx{pass} = flipud(X{pass});
        gy{pass} = flipud(Y{pass});
        TRANSLATIONX{pass} = flipud(U{pass});
        TRANSLATIONY{pass} = flipud(V{pass});
        ROTATION{pass} = flipud(R{pass});
        SCALING{pass} = flipud(S{pass});
        uVal{pass} = flipud(UVAL{pass});
        vVal{pass} = flipud(VVAL{pass});
    end
    
else
    % Initialize the counter for
    % the number of user-specified passes
    piv_pass_number = 1;
    p = 0;
end

% Overwrite the previous pass jobfile and filename info.
FilePaths = FILEPATHS;
JobFile = JOBFILE;

% % Initialize the DWO iteration counter (this is for DWO convergence)
% num_dwo_iterations = 0;

% Loop over the passes.
% for p = 1 : numberOfPasses;
% thisPass is supposed to be the overall pass number
% This is confusing... change to something like intra-pass iteration
while piv_pass_number <= number_of_passes;
    
    % Set this variable which gets saved in the output
%     PASSNUMBER = piv_pass_number;
    
    % Increment the pass counter
    p = p + 1;
    
    % Read deformation flag.
    doImageDeformation = JobFile.Parameters. ...
        Processing(p).Deform.DoImageDeformation;
    
    % Read validation flag
    doValidation = JobFile.Parameters.Processing(p).Validation.DoValidation;
    
    % Read smoothing flag
    doSmoothing = JobFile.Parameters. ...
        Processing(p).Smoothing.DoSmoothing;
    
    % Flag for zero-meaning the interrogation regions
    do_zero_mean = JobFile.Parameters.Processing(p).InterrogationRegion.ZeroMeanRegion;
    
    
    %%% Here we check whether to do image deformation
    % p is the pass number
    % Check deformation flag
    if p > 1 && doImageDeformation
        
        % Smooth field if specified
        if doSmoothing
            
            % Smoothing kernel diameter
            smoothingKernelDiameter = ...
                JobFile.Parameters.Processing(p). ...
                Smoothing.KernelDiameter;
   
            % Smoothing kernel gaussian standard deviation  
            smoothingGaussianStdDev = ...
                JobFile.Parameters.Processing(p). ...
                Smoothing.KernelGaussianStdDev;
    
            % Smooth the velocity field.
            % Not sure how to deal with this 
            % growing field size warning in this case.
            source_field_u{p-1} = smoothField(uVal{p-1}, smoothingKernelDiameter, smoothingGaussianStdDev);
            source_field_v{p-1} = smoothField(vVal{p-1}, smoothingKernelDiameter, smoothingGaussianStdDev);
       
        else
            source_field_u{p-1} = uVal{p-1};
            source_field_v{p-1} = vVal{p-1};
        end
        
        % Create the pixel coordinates.
        [xi_integer, yi_integer] = meshgrid(1:imageWidth, 1:imageHeight);
        
        % Shift the pixel coordinates by 0.5 pixels
        XI = xi_integer - 0.5;
        YI = yi_integer - 0.5;
        
        % Create interpolation structures for the velocity field.
        % Temporary: change from spline to cubic interpolation, and change
        % from lienar to nearest neighbor extrapolation. This is for
        % comparison with prana.
        interpolant_tx = griddedInterpolant(gy{p-1}, gx{p-1}, source_field_u{p-1}, 'cubic', 'nearest');
        interpolant_ty = griddedInterpolant(gy{p-1}, gx{p-1}, source_field_v{p-1}, 'cubic', 'nearest');
        
        % This is the velocity field upsampled to every pixel.
        UI = interpolant_tx(YI, XI);
        VI = interpolant_ty(YI, XI);
        
        % These are the coordinates at which to resample image 1.
        XD1 = XI - UI/2;
        YD1 = YI - VI/2;
        
        % These are the coordinates at which to resample image 2.
        XD2 = XI + UI/2;
        YD2 = YI + VI/2;

        % Resample the images
        image1 = sincBlackmanInterp2(image1_raw, XD1 + 0.5, YD1 + 0.5, 8, 'blackman');
        image2 = sincBlackmanInterp2(image2_raw, XD2 + 0.5, YD2 + 0.5, 8, 'blackman');
        
    else
        
        % If not deform or if we're on the first pass, use the raw images.
        image1 = image1_raw;
        image2 = image2_raw;
    end
    
    % Interrogation region dimensions
    % Region height (pixels)
    regionHeight = JobFile.Parameters. ...
        Processing(p).InterrogationRegion.Height;
    
    % Region width (pixels)
    regionWidth  = JobFile.Parameters. ...
        Processing(p).InterrogationRegion.Width;    
    
    % Fmc difference method
    fmcDifferenceMethod_string = JobFile.Parameters. ...
        Processing(p).FMC.FmcDifferenceMethod; 
    
    % Determine what FMC difference method to use. 
    % fmcDifferenceMethod = 1 is central difference.
    % fmcDifferenceMethod = 2 is forward difference.
    % fmcDifferenceMethod = 3 is backward difference.
    if ~isempty(regexpi(fmcDifferenceMethod_string, 'bac'))
        fmcDifferenceMethod = 3;
    elseif ~isempty(regexpi(fmcDifferenceMethod_string, 'for'))
        fmcDifferenceMethod = 2;
    else
        % Default to central difference
        fmcDifferenceMethod = 1;
    end
    
    % Apodization window parameters.
    spatialWindowFraction = JobFile.Parameters. ...
        Processing(p).InterrogationRegion.SpatialWindowFraction;
    
    % FMC parameters
    fmiWindowSize = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.FMIWindowSize;
    fmiWindowType = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.FMIWindowType;

    % Log-polar Image resampling parameters
    numberOfRings = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.NumberOfRings;
    numberOfWedges = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.NumberOfWedges;
    
    % FFT Parameters
    fftSize = JobFile.Parameters.Processing(p).FFTSize;
    spectrum_height = fftSize(1);
    spectrum_width  = fftSize(2);
    
    rMin = JobFile.Parameters.Processing(p).Correlation.FMC.MinimumRadius;
    rMax = min(spectrum_height, spectrum_width) / 2 - 1;

    % Correlation parameters
    correlationMethod = JobFile.Parameters.Processing(p).Correlation.Method;
    
    % Determine which correlation type to use. Fmc or RPC
    isFmc = ~isempty(regexpi(correlationMethod, 'fmc'));
    isRpc = ~isempty(regexpi(correlationMethod, 'rpc'));
    isScc = ~isempty(regexpi(correlationMethod, 'scc'));
    
    % RPC diameters
    spatialRPCDiameter = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.SpatialRPCDiameter;
    fmiRpcDiameter = JobFile.Parameters. ...
        Processing(p).Correlation.FMC.FMCDiameter; 

    % Create the gaussian intensity window to be applied to the the raw image interrogation regions
    spatialWindow = gaussianWindowFilter_prana(...
        [regionHeight, regionWidth], ...
        spatialWindowFraction .* [regionHeight, regionWidth]);
    
    % Determine the FMI window type.
    isHann1 = ~isempty(regexpi(fmiWindowType, 'hann1'));
    isHann2 = ~isempty(regexpi(fmiWindowType, 'hann2'));
    isGaussianSkew = ~isempty(regexpi(fmiWindowType, 'gauss_skew'));
    
    % Extract the string specifying the subpixel peak fit method
    subpixel_peak_fit_method =...
        JobFile.Parameters.Processing(p).Correlation.PeakFitMethod;
    
    % Peak fit method
    % Convert the subpixel peak-fit method string extracted from the jobfile
    % into the numerical peak-fit method identifier expected by the function
    % subpixel.m
    if ~isempty(regexpi(subpixel_peak_fit_method, 'squ'));
        % Least squares method
         subpixel_peak_fit_method_numerical = 3;
         % Least squares method
    elseif ~isempty(regexpi(subpixel_peak_fit_method, '3'));
        
        % Three-point Gaussian
         subpixel_peak_fit_method_numerical = 1;
    else
        % Three-point Gaussian
         subpixel_peak_fit_method_numerical = 1;
    end
    
    % Create the FMI Window
    if isHann1
        fmiWindow1D = hann1(numberOfRings, [fmiWindowSize(1) fmiWindowSize(2)], fmiWindowSize(3));
        fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
    elseif isHann2
        fmiWindow = hann2([numberOfWedges, numberOfRings], fmiWindowSize(1));
    elseif isGaussianSkew
        fmiWindow1D = gaussianWindowFilter_asymmetric(numberOfRings, fmiWindowFraction);
        fmiWindow = repmat(fmiWindow1D, numberOfWedges, 1);
    else
        fmiWindow = gaussianWindowFilter([numberOfWedges, numberOfRings], fmiWindowSize, 'fraction');
    end

    % Create the gaussian spectral energy filter be applied to the raw image correlation
    imageSpectralFilter = spectralEnergyFilter(regionHeight, regionWidth, spatialRPCDiameter); 

    % Create the FMI spectral filter (i.e. the FMI RPC filter).
    fmiSpectralFilter = spectralEnergyFilter(numberOfWedges, numberOfRings, fmiRpcDiameter); 

    % Make a matrix of the subregion coordinates.
    % Do this only once to increase speed (meshgrid is slow).
    [xImage, yImage] = meshgrid(1 : regionWidth, 1 : regionHeight);

    % Figure out the log polar resampling coordinates.
    % This is done here to avoid calling "meshgrid" 
    % for each interrogation region pair, which sould be slow!
    [xLP, yLP] = LogPolarCoordinates([spectrum_height, spectrum_width], numberOfWedges, numberOfRings, rMin, rMax, 2 * pi);

    % Save the image size to the processing field for easy passing around.
    JobFile.Parameters.Processing(p).Images.Height = imageHeight;
    JobFile.Parameters.Processing(p).Images.Width  = imageWidth;

    % Flag specifying whether or not to do DWO
    doDiscreteWindowOffset = JobFile.Parameters.Processing(p).DWO.DoDiscreteWindowOffset;
    
    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(p).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(p).Grid.Spacing.Y;

    % Make sure the grid buffer is at least half the size of the interrogation region
    gridBufferY = JobFile.Parameters.Processing(p).Grid.Buffer.Y;
    gridBufferX = JobFile.Parameters.Processing(p).Grid.Buffer.X;

    % If the pass number is greater than one, i.e., if at least one pass has
    % finished, and also if discrete window offset is enabled,
    % then interpolate the velocity field from the previous pass
    % onto the grid for the current pass.
    % Round the grid shift values so that grid points are shifted
    % from integer coordinates to integercoordinates
    if p > 1 && doDiscreteWindowOffset
        [gx{p}, gy{p}, gx_01, gy_01, gx_02, gy_02] = ...
            discreteWindowOffset(...
            gx{p-1}, gy{p-1}, ...
            uVal{p-1}, vVal{p-1}, ...
            JobFile.Parameters.Processing(p));
        
    else
         % Generate the list of coordinates that specifies the (X, Y) centers of all of the interrogation regions 
        [ gx{p}, gy{p} ] = gridImage([imageHeight, imageWidth], [gridSpacingY gridSpacingX], gridBufferY, gridBufferX);
   
        % If this is the first pass or if DWO is not specified, then keep
        % the original grid points.
        gx_01 = gx{p};
        gy_01 = gy{p};
        gx_02 = gx{p};
        gy_02 = gy{p}; 
    end
    
    % Calculate the grid shifts for both images. 
    % These should all be integers!!
    gridShiftX_01 = gx_01 - gx{p};
    gridShiftY_01 = gy_01 - gy{p};
    gridShiftX_02 = gx_02 - gx{p};
    gridShiftY_02 = gy_02 - gy{p};

    % Determine the size of the grid (number of grid points)..
    [num_grid_rows, num_grid_cols] = size(gx_01);

    % Determine the number of interrogation regions to be correlated
    num_regions = num_grid_rows * num_grid_cols;

    % Extract the subregions from the first image.
    regionMatrix1 = extractSubRegions(...
        image1, [regionHeight, regionWidth], ...
        gx_01(:), gy_01(:));
    
    % Extract the subregions from the second image.
    regionMatrix2 = extractSubRegions(...
        image2, [regionHeight, regionWidth], ...
        gx_02(:), gy_02(:));
    
    % Preallocate memory for the vectors to hold the estimates of translation, rotation, and scaling.
    estimatedTranslationY = zeros(num_regions, 1);
    estimatedTranslationX = zeros(num_regions, 1);
    estimatedRotation = zeros(num_regions, 1); 
    estimatedScaling = ones(num_regions, 1);

    % Initialize RPC peak height ratio vector.
    spatialPeakRatio = zeros(num_regions, 1);
    
    % Start a timer
    t = tic;
    
    % Do all the correlations for the image.
    for k = 1 : num_regions
        
        fprintf(1, [num2str(k) ' of ' num2str(num_regions) '\n']);
        
        % Extract the subregions from the subregion stacks.
        subRegion1 = regionMatrix1(:, :, k);
        subRegion2 = regionMatrix2(:, :, k);
        
        % Zero mean the region if requested
        if do_zero_mean
            subRegion1 = zero_mean_region(subRegion1);
            subRegion2 = zero_mean_region(subRegion2);
        end
                
        % Perform FMC processing. 
        if isFmc
            % Perform the FMC correlation.
            [estimatedTranslationY(k), estimatedTranslationX(k),...
            estimatedRotation(k), estimatedScaling(k), ...
            ~, spatialPeakRatio(k), ...
            ~, ~] = ...
            ...
            FMC(subRegion1, subRegion2, spatialWindow, imageSpectralFilter,...
            fmiWindow, fmiSpectralFilter, ...
            xImage, yImage, ...
            spectrum_height, spectrum_width,...
            xLP, yLP, rMin, rMax, fmcDifferenceMethod, COMPILED);
        
        % Perform RPC analysis 
        % The zero in this input means "Do not search multiple peaks,"
        % i.e., use only the primary peak.
        elseif isRpc
            [estimatedTranslationY(k), estimatedTranslationX(k), ...
                rpcPlane, ~, ~]...
                = RPC(...
                spatialWindow .* subRegion1, ...
                spatialWindow .* subRegion2,...
                imageSpectralFilter, subpixel_peak_fit_method_numerical); 

            % Measure the peak height ratio
            spatialPeakRatio(k) = measurePeakHeightRatio(rpcPlane, COMPILED);
          
        % Perform SCC analysis.
        elseif isScc
            [estimatedTranslationY(k), estimatedTranslationX(k), ...
                scc_plane]...
                = SCC(...
                spatialWindow .* subRegion1, ...
                spatialWindow .* subRegion2,...
                subpixel_peak_fit_method_numerical);
            
                % Measure the peak height ratio
                spatialPeakRatio(k) = measurePeakHeightRatio(scc_plane, COMPILED);
        end
             
    end % end for k = 1 : nRegions
        
    % Inform the user
    disp(['Correlation times (pass ' num2str(p) '): ' num2str(toc(t)) ' sec' ]);
    disp('');
    
    % Reshape the raw measured displacements into matrices.
    tx_raw{p} = reshape(estimatedTranslationX, num_grid_rows, num_grid_cols);
    ty_raw{p} = reshape(estimatedTranslationY, num_grid_rows, num_grid_cols);
    
    % Shift the measured velocities by the deform or DWO values. 
    if doImageDeformation
        
        % Only shift for the second and greater iterations
        if p > 1
            
            % Temporary change to cubic and nearest
            % Create interpolant structures
              interpolant_tx = griddedInterpolant(YI, XI, UI, 'cubic', 'nearest');
              interpolant_ty = griddedInterpolant(YI, XI, VI, 'cubic', 'nearest');

            % Evaluate the interpolant structures
              tx_shift = interpolant_tx(gy{p}, gx{p});
              ty_shift = interpolant_ty(gy{p}, gx{p});

        else
             % For the first pass, set the shift values to zero.
              tx_shift = zeros(size(gx{p}));
              ty_shift = zeros(size(gx{p}));
        end
        
        % Add the shift values to the measured displacements for deform.
        TRANSLATIONX{p} = tx_raw{p} + tx_shift;
        TRANSLATIONY{p} = ty_raw{p} + ty_shift;
        
    else     
        % If not deform then just shift the velocities by the DWO grid
        % shift values, which are already zeros if DWO wasn't used.
        TRANSLATIONX{p} = tx_raw{p} + gridShiftX_02 - gridShiftX_01;
        TRANSLATIONY{p} = ty_raw{p} + gridShiftY_02 - gridShiftY_01;
    end
    
    % Reshape the rotation and scaling measurements into matrices.
    ROTATION{p} = reshape(estimatedRotation, num_grid_rows, num_grid_cols);
    SCALING{p} =  reshape(estimatedScaling, num_grid_rows, num_grid_cols);

    % Reshape the peak ratio measurements into matrices.
    SPATIAL_PEAK_RATIO{p} = flipud(reshape(spatialPeakRatio, num_grid_rows, num_grid_cols));

    % Run Prana's validation code. Note that right now the rotation
    % The previous codes to do this were "universalOutlierDetection.m"
    % and "universalOutlierReplacement.m"
    % estimate isn't validated.
    if doValidation
        
        % Validation parameters
        val_parameters = JobFile.Parameters.Processing(p).Validation;
        
        % Perform the validation
        [uVal{p}, vVal{p}, isOutlier{p}] = validateField_prana(...
        gx{p}, gy{p}, ...
        TRANSLATIONX{p}, TRANSLATIONY{p}, ...
        val_parameters);
    else
        uVal{p} = zeros(size(tx_raw{p}));
        vVal{p} = zeros(size(ty_raw{p}));
        isOutlier{p} = zeros(size(tx_raw{p}));     
    end
    
    % Read the flag for convergence checking
    check_for_convergence = ...
        JobFile.Parameters.Processing(p).CheckConvergence;
    
    % Check convergence?
    % Dummy variables for now
    if check_for_convergence && p > 1
        has_converged(p) = check_convergence(u{p}, v{p}, ...
            u{p-1}, v{p-1}, ...
            convergence_criterion);
        
        if has_converged(p) || num_iterations > max_iterations
            
            % Increment the pass number.
            piv_pass_number = piv_pass_number + 1;
            
            % Reset the number of iterations
            num_iterations = 0;
        else
            
            % Increment the number of iterations performed
            % during this pass
            num_iterations = num_iterations + 1;
            
            % Update the jobfile
            JobFile.Parameters.Processing(p + 1 : end + 1) = ...
                JobFile.Parameters.Processing(p : end);
        end
        
        piv_pass_number = piv_pass_number + 1;
        has_converged(p) = 0;
        
    end
    
        
end % End for p = 1 : numberOfPasses  

% Number of passes that ended up getting run
finalNumberOfPasses = length(JobFile.Parameters.Processing);

% This flips everything over the Y axis and saves the variables to the
% output structures.
for p = 1 : finalNumberOfPasses
    X{p} = flipud(gx{p});
    Y{p} = flipud(gy{p});
    U{p} = flipud(TRANSLATIONX{p});
    V{p} = flipud(TRANSLATIONY{p});
    R{p} = flipud(ROTATION{p});
    S{p} = flipud(SCALING{p});
    UVAL{p} = flipud(uVal{p});
    VVAL{p} = flipud(vVal{p});
    IS_OUTLIER{p} = flipud(isOutlier{p});
    
    % TEMPORARY: Don't validate rotation estimate.
    RVAL{p} = zeros(size(R{p}));
    
end

% Rename variable
CONVERGED = hasConverged;

% Set source field variables to zeros if only one pass was specifed.
if number_of_passes < 2
    source_field_u{1} = zeros(size(X{1}));
    source_field_v{1} = zeros(size(X{1}));
end

% Save the results
save(FilePaths.OutputFilePath, ...
    'X', 'Y', 'U', 'V', 'R', 'S', 'IS_OUTLIER',...
    'UVAL', 'VVAL', 'RVAL', 'tx_raw', 'ty_raw',...
    'SPATIAL_PEAK_RATIO', 'PASSNUMBER', 'CONVERGED', ...
    'FilePaths', 'JobFile');

end






















