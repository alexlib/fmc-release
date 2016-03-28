
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

% Chop the jobfile processing parameters down to the first n passes
JOBFILE.Parameters.Processing = JOBFILE.Parameters.Processing(1 : number_of_passes);

% Read in the clean images.
image1_import = double(imread(FILEPATHS.FirstImagePath));
image2_import = double(imread(FILEPATHS.SecondImagePath));

% Image dimensions
[imageHeight, imageWidth, num_channels] = size(image1_import);

% Save the image size to the processing field for easy passing around.
for p = 1 : number_of_passes
    JOBFILE.Parameters.Processing(p).Images.Height = imageHeight;
    JOBFILE.Parameters.Processing(p).Images.Width  = imageWidth;
end

% Loads the mask or create a mask of ones if no mask is specified.
do_mask = JOBFILE.Parameters.Mask.DoMasking;

% Default to a mask of ones (i.e., unmasked)
mask = ones(imageHeight, imageWidth);

% Load the mask if requested
if do_mask
    mask_path = JOBFILE.Parameters.Mask.Path;
    if exist(mask_path, 'file')
        mask = double(imread(mask_path));
    else
        fprintf('Mask not found!\n');
    end
end
    
% Extract the appropriate color channel. channel
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
    
else
    % Default to taking the first channel if none is specified.
    image1_raw = image1_import(:, :, 1);
    image2_raw = image2_import(:, :, 1);
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
        thisPass = 2;
    else
        thisPass = JOBFILE.JobOptions.StartPass;
    end
    
    % This sets a loop counter for the PIV passes.
    p = thisPass - 1;
    
    % Create the function variables from the saved data
    % This just entails flipping all the coordinate and 
    % displacement data
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
    % Initialize the counter for the number of user-specified passes
    thisPass = 1;
    p = 0;
end

% Overwrite the previous pass jobfile and filename info.
FilePaths = FILEPATHS;
JobFile = JOBFILE;

% Initialize the iteration counter 
% (this is for DWO and deform convergence)
num_iterations = 0;

% Loop over the passes.
% for p = 1 : numberOfPasses;
while thisPass <= number_of_passes;
    %% This region of the code reads the parameters for each pass.
    
    % Increment the pass counter 
    p = p + 1;
    
    % Write the number of the specified pass that's currently executing
    PASSNUMBER(p) = thisPass;
    
    % Read smoothing parameters
    doSmoothing = JobFile.Parameters.Processing(p).Smoothing.DoSmoothing;
    
    % Smoothing kernel diameter
    smoothingKernelDiameter = JobFile.Parameters.Processing(p).Smoothing.KernelDiameter;
   
    % Smoothing kernel gaussian standard deviation  
    smoothingGaussianStdDev = JobFile.Parameters.Processing(p).Smoothing.KernelGaussianStdDev;
    
    % Flag for zero-meaning the interrogation regions
    do_zero_mean = JobFile.Parameters.Processing(p).InterrogationRegion.ZeroMeanRegion;
    
    % Interrogation region dimensions
    regionHeight = JobFile.Parameters.Processing(p).InterrogationRegion.Height;
    regionWidth  = JobFile.Parameters.Processing(p).InterrogationRegion.Width;    
        
    % FFT Parameters
    fftSize = JobFile.Parameters.Processing(p).FFTSize;
    spectrum_height = fftSize(1);
    spectrum_width  = fftSize(2);
    
    %% These things are specific to FMC, but run for all methods
    % because of some issues with the for loop (I think).
    
    % Window sizes and types for FMC; this is specific to FMC.
    fmiWindowSize = ...
        JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowSize;
    fmiWindowType = ...
        JobFile.Parameters.Processing(p).InterrogationRegion.FMIWindowType;

    % Fmc difference method
    fmcDifferenceMethod_string =...
        JobFile.Parameters.Processing(p).FMC.FmcDifferenceMethod; 
    
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

    % Image resampling parameters
    % This is specific to FMC processing.
    numberOfRings = JobFile.Parameters.Processing(p).Correlation.FMC.NumberOfRings;
    numberOfWedges = JobFile.Parameters.Processing(p).Correlation.FMC.NumberOfWedges;
    rMin = JobFile.Parameters.Processing(p).Correlation.FMC.MinimumRadius;
    rMax = min(spectrum_height, spectrum_width) / 2 - 1;
    
    % Figure out the log polar resampling coordinates.
    % This is specific to FMC.
    % It's dumb to do this even if FMC isn't specified, but
    % I haven't taken the time to fix it. 
    [xLP, yLP] = LogPolarCoordinates([spectrum_height, spectrum_width], ...
        numberOfWedges, numberOfRings, rMin, rMax, 2 * pi);
    
    % "RPC" diameter used for the FMC correlations.
    fmiRpcDiameter = ...
        JobFile.Parameters.Processing(p).Correlation.FMC.FilterDiameter; 
    
    % Create the FMI spectral filter (i.e. the FMI RPC filter).
    fmiSpectralFilter = spectralEnergyFilter(...
        numberOfWedges, numberOfRings, fmiRpcDiameter); 

    % Determine the FMI window type. This is specific to FMC.
    isHann1 = ~isempty(regexpi(fmiWindowType, 'hann1'));
    isHann2 = ~isempty(regexpi(fmiWindowType, 'hann2'));
    isGaussianSkew = ~isempty(regexpi(fmiWindowType, 'gauss_skew'));
    
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
    
    %% This section reads and specifies some more parameters
    % that aren't specific to FMC.
    
    % RPC diameters
    spatialRPCDiameter = JobFile.Parameters.Processing(p).Correlation.RPC.FilterDiameter;
    
    % Apodization window parameters (window size for the interrogation
    % regions; applies to all methods).
    spatialWindowFraction = ...
        JobFile.Parameters.Processing(p).InterrogationRegion.SpatialWindowFraction;
      
    % Create the gaussian intensity window to be applied
    % to the raw image interrogation regions
    spatialWindow = gaussianWindowFilter_prana(...
        [regionHeight, regionWidth], ...
        spatialWindowFraction .* [regionHeight, regionWidth]);
    
    % Extract the string specifying the subpixel peak fit method
    subpixel_peak_fit_method =...
        JobFile.Parameters.Processing(p).Correlation.PeakFitMethod;
    
    % Peak fit method
    % Convert the subpixel peak-fit method string extracted from the jobfile
    % into the numerical peak-fit method identifier expected by the function
    % subpixel.m
    if ~isempty(regexpi(subpixel_peak_fit_method, 'squ'));
         subpixel_peak_fit_method_numerical = 3;
    elseif ~isempty(regexpi(subpixel_peak_fit_method, '3'));
         subpixel_peak_fit_method_numerical = 1;
    else
         subpixel_peak_fit_method_numerical = 1;
    end

    % This creates the RPC filter that gets applied to the correlations
    % of the particle images; it applies to RPC and FMC (and SPC later).
    imageSpectralFilter = spectralEnergyFilter(...
        regionHeight, regionWidth, spatialRPCDiameter); 

    % Make a matrix of the subregion coordinates.
    % Do this only once to increase speed (meshgrid is slow).
    [xImage, yImage] = meshgrid(1 : regionWidth, 1 : regionHeight);
    
    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(p).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(p).Grid.Spacing.Y;

    % Make sure the grid buffer is at least half the size of the interrogation region
    gridBufferY = JobFile.Parameters.Processing(p).Grid.Buffer.Y;
    gridBufferX = JobFile.Parameters.Processing(p).Grid.Buffer.X;
    
    %% This section handles modifying the images and creating or
    % modifying the grid based on previous iterations if those options 
    % are specified (e.g., image deformation, discrete window offset)
    
    % Iterative method (DWO, Deform, etc)
    iterative_method = JobFile.Parameters.Processing(p).Iterative.Method;
    
    % This determines whether DWO convergence iterations were specified
    converge_iterative_method = JobFile.Parameters.Processing(p).Iterative.Converge;
    
    % Check what iterative method is requested (if any)
    if regexpi(iterative_method, 'def') % This is the case for deform.
        doImageDeformation = true;
        doDiscreteWindowOffset = false;
    elseif regexpi(iterative_method, 'dwo') % Case for discrete window offset
        doImageDeformation = false;
        doDiscreteWindowOffset = true;
    else % Case for neither DWO nor Deform
        doImageDeformation = false;
        doDiscreteWindowOffset = false;
        
        % Don't try to converge anything if
        % no iterative method was specified.
        converge_iterative_method = false;
    end
    
    % These lines determine what velocity field to use
    % as the source-field for deform, DWO, etc.
    % Smoothing happens here. An important
    % implication of smoothing happening here (and
    % not saving the smoothed fields as their own variables)
    % is that smoothed data are not saved to disk. 
    % This is so that the user never gets a field that's unintentinoally
    % smoothed; it kind of forces you to look at the unsmoothed data.
    % Smoothing is easy in post processing. However, the smoothing
    % still takes place in the code so that deform, and DWO are stable.
    if p > 1
        % Smooth field if specified
        if doSmoothing
            % Smooth the velocity field.
            source_field_u{p-1} = smoothField(uVal{p-1}, ...
                smoothingKernelDiameter, smoothingGaussianStdDev);
            source_field_v{p-1} = smoothField(vVal{p-1},...
                smoothingKernelDiameter, smoothingGaussianStdDev);
       
        else
            source_field_u{p-1} = uVal{p-1};
            source_field_v{p-1} = vVal{p-1};
        end
    end
    
    % Perform image deformation if requested.
    % The p > 1 statement says to only deform the images
    % (i.e., don't try to smooth before any
    % velocities have been calculated).
    if p > 1 && doImageDeformation
        
        % Start a timer to measure the number of 
        % seconds spent on image deformation.
        deform_tic = tic;
 
        % Create the pixel coordinates.
        [xi_integer, yi_integer] = meshgrid(1:imageWidth, 1:imageHeight);
        
        % Shift the pixel coordinates by 0.5 pixels
        XI = xi_integer - 0.5;
        YI = yi_integer - 0.5;
        
        % Create interpolation structures for the velocity field.
        % Temporary: change from spline to cubic interpolation, and change
        % from lienar to nearest neighbor extrapolation. This is for
        % comparison with prana.
        interpolant_tx = scatteredInterpolant(gy{p-1}, gx{p-1}, ...
            source_field_u{p-1}, 'nearest', 'nearest');
        interpolant_ty = scatteredInterpolant(gy{p-1}, gx{p-1}, ...
            source_field_v{p-1}, 'nearest', 'nearest');

        % This is the velocity field upsampled to every pixel.
        UI = interpolant_tx(YI, XI);
        VI = interpolant_ty(YI, XI);
        
        % These are the coordinates at which to resample image 1.
        XD1 = XI - UI/2;
        YD1 = YI - VI/2;
        
        % These are the coordinates at which to resample image 2.
        XD2 = XI + UI/2;
        YD2 = YI + VI/2;

        % Deform (resample) the first image
        fprintf('Deforming image 1...\n')
        image1 = sincBlackmanInterp2(image1_raw, XD1 + 0.5, YD1 + 0.5, 8, 'blackman');
        
        % Deform (resample) the first image
        fprintf('Deforming image 2...\n')
        image2 = sincBlackmanInterp2(image2_raw, XD2 + 0.5, YD2 + 0.5, 8, 'blackman');
        
        % End the deform timer.
        deform_toc = toc(deform_tic);
        
        % Display the number of seconds spent on image deformation.
        fprintf(1, 'Deform time: %0.2f seconds.\n', deform_toc);
        
    else
        % If not deform or if we're on the first pass, use the raw images.
        image1 = image1_raw;
        image2 = image2_raw;
    end
      
    % If the pass number is greater than one, i.e., if at least one pass has
    % finished, and also if discrete window offset is enabled,
    % then interpolate the velocity field from the previous pass
    % onto the grid for the current pass.
    % Round the grid shift values so that grid points are shifted
    % from integer coordinates to integercoordinates
    if p > 1 && doDiscreteWindowOffset
        [gx{p}, gy{p}, gx_01, gy_01, gx_02, gy_02] = ...
            discreteWindowOffset(gx{p-1}, gy{p-1}, ...
            source_field_u{p-1}, source_field_v{p-1}, ...
            JobFile.Parameters.Processing(p));
        
    else
         % Generate the list of coordinates that specifies the (X, Y) 
         % centers of all of the interrogation regions 
        [ gx{p}, gy{p} ] = gridImage([imageHeight, imageWidth], ...
            [gridSpacingY gridSpacingX], ...
            gridBufferY, gridBufferX);
   
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
    [numRows, numColumns] = size(gx_01);

    % Determine the number of interrogation regions to be correlated
    nRegions = numRows * numColumns;
    
    % Extract the subregions from image 1.
    regionMatrix1 = extractSubRegions(image1 .* mask, ...
        [regionHeight, regionWidth], gx_01(:), gy_01(:));
    
    % Extract the subregions from image 2.
    regionMatrix2 = extractSubRegions(image2 .* mask, ...
        [regionHeight, regionWidth], gx_02(:), gy_02(:));
    
    % Preallocate memory for the vectors to hold the 
    % estimates of translation, rotation, and scaling.
    % Some of these don't get saved right now; they're here
    % in case we decide to change this.
    estimatedTranslationY = zeros(nRegions, 1);
    estimatedTranslationX = zeros(nRegions, 1);
    estimatedRotation = zeros(nRegions, 1); 
    estimatedScaling = ones(nRegions, 1);
   
    % Initialize FMC peak height ratio vector.
    fmcPeakRatio = zeros(nRegions, 1);

    % Initialize RPC peak height ratio vector.
    spatialPeakRatio = zeros(nRegions, 1);
    
    % Start a timer
    t = tic;
    
    %% This section does the correlations for the current and iteration.
    
    % Read the correlation method
    correlationMethod = ...
        JobFile.Parameters.Processing(p).Correlation.Method;
    
    % Determine which correlation type to use (FMC, RPC, SCC)
    isFmc = ~isempty(regexpi(correlationMethod, 'fmc'));
    isRpc = ~isempty(regexpi(correlationMethod, 'rpc'));
    isScc = ~isempty(regexpi(correlationMethod, 'scc'));
    
    % Do all the correlations for the image.
    parfor k = 1 : nRegions
        
        % This line prints to screen for every single
        % correlation; it's useful if correlations are 
        % intensive and taking a long time, and you want 
        % to keep tabs on what's happening rather than 
        % just staring at a blank screen. Uncomment it if you
        % want a more "verbose" output.
%         fprintf(1, 'On region %d of %d\n', k, nRegions);
        
        % Extract the subregions from the subregion stacks.
        subRegion1 = regionMatrix1(:, :, k);
        subRegion2 = regionMatrix2(:, :, k);
        
        if max(subRegion1(:)) > 0
        
            % Zero mean the region if requested
            if do_zero_mean
                subRegion1 = zero_mean_region(subRegion1);
                subRegion2 = zero_mean_region(subRegion2);
            end

            % Perform FMC correlation. 
            if isFmc
                % Perform the FMC correlation.
                [estimatedTranslationY(k), estimatedTranslationX(k),...
                estimatedRotation(k), estimatedScaling(k), ...
                fmcPeakRatio(k), spatialPeakRatio(k)] = ...
                ...
                FMC(subRegion1, subRegion2, spatialWindow, imageSpectralFilter,...
                fmiWindow, fmiSpectralFilter, ...
                xImage, yImage, ...
                spectrum_height, spectrum_width,...
                xLP, yLP, rMin, rMax, fmcDifferenceMethod, COMPILED);

            % Perform RPC correlation 
            % The zero in this input means "Do not search multiple peaks,"
            % i.e., use only the primary peak.
            elseif isRpc
                [estimatedTranslationY(k), estimatedTranslationX(k), ...
                    rpcPlane]...
                    = RPC(...
                    spatialWindow .* subRegion1, ...
                    spatialWindow .* subRegion2,...
                    imageSpectralFilter, ...
                    subpixel_peak_fit_method_numerical, COMPILED); 

                % Measure the peak height ratio
                spatialPeakRatio(k) = measurePeakHeightRatio(rpcPlane, COMPILED);

            % Perform SCC analysis.
            elseif isScc
                [estimatedTranslationY(k), estimatedTranslationX(k)]...
                    = SCC(...
                    spatialWindow .* subRegion1, ...
                    spatialWindow .* subRegion2,...
                    subpixel_peak_fit_method_numerical, COMPILED);
            end
        end
        
    end % end for k = 1 : nRegions
        
    % Inform the user
    fprintf('Correlation time: %0.2f seconds.\n', toc(t));
    
    % Reshape the raw measured displacements into matrices.
    tx_raw = reshape(estimatedTranslationX, numRows, numColumns);
    ty_raw = reshape(estimatedTranslationY, numRows, numColumns);
    
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
        TRANSLATIONX{p} = tx_raw + tx_shift;
        TRANSLATIONY{p} = ty_raw + ty_shift;
        
    else     
        % If not deform then just shift the velocities by the DWO grid
        % shift values, which are already zeros if DWO wasn't used.
        TRANSLATIONX{p} = tx_raw + gridShiftX_02 - gridShiftX_01;
        TRANSLATIONY{p} = ty_raw + gridShiftY_02 - gridShiftY_01;
    end
    
    % Reshape the rotation and scaling measurements into matrices.
    ROTATION{p} = reshape(estimatedRotation, numRows, numColumns);
    SCALING{p} =  reshape(estimatedScaling, numRows, numColumns);

    % Reshape the peak ratio measurements into matrices.
    SPATIAL_PEAK_RATIO{p} = flipud(reshape(spatialPeakRatio, numRows, numColumns));
    FMC_PEAK_RATIO{p} = flipud(reshape(fmcPeakRatio, numRows, numColumns));
    
    % Extract the validation parameters
    validation_parameters = JobFile.Parameters.Processing(p).Validation;
    
    % Run Prana's validation code. Note that right now the rotation
    % The previous codes to do this were "universalOutlierDetection.m"
    % and "universalOutlierReplacement.m"
    % estimate isn't validated.
    [uVal{p}, vVal{p}, isOutlier{p}] = ...
        validateField_prana(gx{p}, gy{p}, ...
        TRANSLATIONX{p}, TRANSLATIONY{p}, ...
        validation_parameters);
        
    % Check for convergence if it's requested. 
    % If the velocity estimate has converged, go on
    % to the next user-specified pass. Otherwise, 
    % repeat the previous pass. 
    
    % If DWO convergence was specified (and at least one pass has
    % completed), then check the other parameters regarding convergence
    if converge_iterative_method
        % Inform the user
        disp('Convergence requested. Checking convergence.')
        
        % Determine the maximum number of iterations specified
        max_iterations = ...
            JobFile.Parameters.Processing(p).Iterative.MaxIterations;
        
        % Determine the iterative method convergence criteria
        convergence_criterion = ...
            JobFile.Parameters.Processing(p).Iterative.ConvergenceCriterion;
        
        % If at least one iteration has been completed...
        if num_iterations > 0
            % Determine the 2-norm of the velocity field components compared to
            % the previous pass. This is the metric against which the
            % convergence criteria is compared.
            
            % Find the indices where all the displacements are nonzero.
            % Identically zero displacements should only evaluate
            % where the image was zero; i.e., masked.
            inds = ...
                abs(TRANSLATIONX{p}(:)) > 0 & ...
                abs(TRANSLATIONY{p}(:)) > 0 & ...
                abs(TRANSLATIONX{p - 1}(:)) > 0 & ...
                abs(TRANSLATIONY{p - 1}(:)) > 0;
            
            % Calculate the norm of the magnitude of
            % the difference of the horizontal velocities. 
            uNorm(p) = ...
                mean(abs(TRANSLATIONX{p}(inds) - TRANSLATIONX{p-1}(inds)));
            
            % Calculate the norm of the magnitude of
            % the difference of the vertical velocities. 
            vNorm(p) = ...
                mean(abs(TRANSLATIONY{p}(inds) - TRANSLATIONY{p-1}(inds)));
            
            % Inform the user
            disp(...
                ['U norm: ' num2str(uNorm(p),...
                '%10.3e') '    V norm: ' ...
                num2str(vNorm(p), '%10.3e') '    Criteria: '...
                num2str(convergence_criterion, '%10.3e')]);

            % Check if the convergence criteria have been reached
            hasConverged(p) = ...
                min(uNorm(p) <= convergence_criterion , ...
                vNorm(p) <= convergence_criterion);
        else
            % Convergence is never reached before the first iteration.
            hasConverged(p) = 0;
        end
        
        % If the velocity field has converged or
        % the max number of iterations has been reached
        if hasConverged(p)
           
            % Inform the user that the pass converged
            fprintf('Pass %d converged after %d iterations of %s.\n',...
                thisPass, num_iterations, iterative_method)
            
            % Save the number of iterations
            iterations(thisPass) = num_iterations;
     
            % Reset the Iteration counter
            num_iterations = 0;
            
            % Increment the counter for the user-specified passes
            thisPass = thisPass + 1;
            
            % Inform the user that the pass will now increment.
            if thisPass < number_of_passes
                fprintf('Incrementing pass.\n\n');
            else
                fprintf('\n');
            end
           
        % Max iterations reached.
        elseif num_iterations > max_iterations
            % Inform the user that the max number of iterations was reached
            fprintf(...
                'Max number of iterations reached for pass %d.\n', thisPass);
            
            % Save the number of iterations
            iterations(thisPass) = num_iterations;
     
            % Reset the Iteration counter
            num_iterations = 0;
            
            % Increment the counter for the user-specified passes
            thisPass = thisPass + 1;
        
        % Neither convergence nor max iterations have been reached.    
        else
            
            % Inform the user
            fprintf(...
                'Pass %d has not converged after %d iterations of %s.\n\n', ...
                thisPass, num_iterations, iterative_method);
            fprintf(1, 'Iterating %s...\n', iterative_method);

            % Increment the iteration counter
            num_iterations = num_iterations + 1; 

            % Update the jobfile with a new pass
            % Then exit the convergence-checking loop with the updated jobfile
            JobFile.Parameters.Processing(p + 1 : end + 1) = JobFile.Parameters.Processing(p : end);
            
        end
        
    else
        % Increment the counter for the user-specified passes
        thisPass = thisPass + 1;
        hasConverged(p) = 0;
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

% Save the results to file.
save(FilePaths.OutputFilePath, ...
    'X', 'Y', 'U', 'V', 'R', 'S', 'IS_OUTLIER',...
    'UVAL', 'VVAL', 'RVAL',...
    'PASSNUMBER', ...
    'hasConverged', ...
    'FilePaths', 'JobFile');

end






















