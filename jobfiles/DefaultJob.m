function defaultJob = DefaultJob();

% Job options
defaultJob.JobOptions.NumberOfPasses = 1;
defaultJob.JobOptions.SkipExisting = 0;
defaultJob.JobOptions.StartFromExistingField = 0;
defaultJob.JobOptions.StartPass = 1;
defaultJob.JobOptions.RunCompiled = false;

% Image parameters
defaultJob.Parameters.Images.Directory = '~/Desktop/piv_images/raw';
defaultJob.Parameters.Images.BaseName = 'lambvortex_h1024_w1024_';
defaultJob.Parameters.Images.Extension = '.tiff';
defaultJob.Parameters.Images.NumberOfDigits = 6;
defaultJob.Parameters.Images.CorrelationStep = 3;

% Start and end images
defaultJob.Parameters.Images.Start = 1;
defaultJob.Parameters.Images.End = 1;
defaultJob.Parameters.Images.FrameStep = 1;
defaultJob.Parameters.Images.ColorChannel = 1;

% Output (vector) parameters
defaultJob.Parameters.Vectors.Directory = '~/Desktop/piv_images/vect';
defaultJob.Parameters.Vectors.BaseName = 'frame_';
defaultJob.Parameters.Vectors.NumberOfDigits = 6;

% Masking
defaultJob.Parameters.Mask.DoMasking = true;
defaultJob.Parameters.Mask.Path = '~/Desktop/mask.tif';

% Grid parameters
% These are the same for all processing methods
defaultJob.Parameters.Processing.Grid.Spacing.X = 32; % Horizontal spacing between grid points (pixels)
defaultJob.Parameters.Processing.Grid.Spacing.Y = 32; % Vertical spacing between grid points (pixels)
defaultJob.Parameters.Processing.Grid.Buffer.Y = [0 0]; % Top and bottom grid buffer (pixels)
defaultJob.Parameters.Processing.Grid.Buffer.X = [0 0]; % Left and right grid buffer (pixels)

% Interrogation Region Parameters
% These are the same for all processing methods
defaultJob.Parameters.Processing.InterrogationRegion.Height = 128;
defaultJob.Parameters.Processing.InterrogationRegion.Width = 128;
defaultJob.Parameters.Processing.InterrogationRegion.SpatialWindowFraction = [0.5 0.5];
defaultJob.Parameters.Processing.InterrogationRegion.ZeroMeanRegion = 1;

% FMC options
defaultJob.Parameters.Processing.InterrogationRegion.FMIWindowSize = [2, 2, 0];
defaultJob.Parameters.Processing.InterrogationRegion.FMIWindowType = 'hann1';

% Correlation parameters
% Correlation type
defaultJob.Parameters.Processing.Correlation.Method = 'rpc'; 

% Peak fit method ('least-squares', ...)
defaultJob.Parameters.Processing.Correlation.PeakFitMethod =...
    '3-point';

% FMC parameters
% FMC Image resampling parameters
defaultJob.Parameters.Processing.Correlation.FMC.NumberOfRings = 64;
defaultJob.Parameters.Processing.Correlation.FMC.NumberOfWedges = 256;
defaultJob.Parameters.Processing.Correlation.FMC.MinimumRadius = 2;

% FMC parameters for windowing the FMI log polar images.
defaultJob.Parameters.Processing.Correlation.FMC.FMIWindowSize = [2, 2, 0];
defaultJob.Parameters.Processing.Correlation.FMC.FMIWindowType = 'hann1';

% RPC spectral filter diameter for the FMC correlation (pixels)
defaultJob.Parameters.Processing.Correlation.FMC.FilterDiameter = 3.3; 
defaultJob.Parameters.Processing.Correlation.FMC.FMCFilterType = 'relative';

% Spatial RPC diameter
defaultJob.Parameters.Processing.Correlation.RPC.FilterDiameter = 2.8;


% FFT Dimensions.
% Is this specific to FMC? 
defaultJob.Parameters.Processing.FFTSize = [128, 128];

% Search multiple peaks?
defaultJob.Parameters.Processing.MultiPeak = false;

% Get rid of this?
defaultJob.Parameters.Processing.FMC.FmcDifferenceMethod = 'forward';

%%%%%%%%% 

% Iterative method paramters (DWO and deform)
% Method can be 'deform', 'dwo', or 'none'
defaultJob.Parameters.Processing.Iterative.Method = 'dwo';

% Attempt to converge the iterative method?
defaultJob.Parameters.Processing.Iterative.Converge = true;

% Maximum number of iterations to use for the pass.
defaultJob.Parameters.Processing.Iterative.MaxIterations = 4;

% Convergence criterion for the pass' iterative method
% (average change in both components of the vector magnitude
% between iterations)
defaultJob.Parameters.Processing.Iterative.ConvergenceCriterion = 0.01;

%%%%%%%%%%%%%

% Smoothing parameters
defaultJob.Parameters.Processing.Smoothing.DoSmoothing = false;
defaultJob.Parameters.Processing.Smoothing.KernelDiameter = 7; 
defaultJob.Parameters.Processing.Smoothing.KernelGaussianStdDev = 1;

% Universal Outlier Detection Parameters
defaultJob.Parameters.Processing.Validation.DoValidation = true;
defaultJob.Parameters.Processing.Validation.UodStencilRadius = 1;
defaultJob.Parameters.Processing.Validation.UodThreshold = 3;
defaultJob.Parameters.Processing.Validation.UodExpectedDifference = [0.1, 0.1];
defaultJob.Parameters.Processing.Validation.UodMedianThreshold = [3, 2];
defaultJob.Parameters.Processing.Validation.UodWindowSize = [3, 3; 3, 3];
defaultJob.Parameters.Processing.Validation.UThresh = [-inf, inf];
defaultJob.Parameters.Processing.Validation.VThresh = [-inf, inf];

% Default Processing parameters
defaultProcessing = defaultJob.Parameters.Processing;

end
