function JOBLIST = fmcJobList()

% Job options
DefaultJob.JobOptions.NumberOfProcessors = 1;
DefaultJob.JobOptions.ImageRotationAngle = 0;
DefaultJob.JobOptions.NumberOfPasses = 1;
DefaultJob.JobOptions.SkipExisting = 0;
DefaultJob.JobOptions.LeftHanded = 0;
% DefaultJob.JobOptions.ComparisonType = 'Eulerian';
DefaultJob.JobOptions.StartFromExistingField = 0;
DefaultJob.JobOptions.StartPass = 1;
DefaultJob.JobOptions.RunCompiled = 1;

% Image parameters
DefaultJob.Parameters.Images.Directory = '~/Desktop/piv_images/raw';
DefaultJob.Parameters.Images.BaseName = 'lambvortex_h1024_w1024_';
DefaultJob.Parameters.Images.Extension = '.tiff';
DefaultJob.Parameters.Images.NumberOfDigits = 6;
DefaultJob.Parameters.Images.CorrelationStep = 1;

% Stard and end images
DefaultJob.Parameters.Images.Start = 1;
DefaultJob.Parameters.Images.End = 1;
DefaultJob.Parameters.Images.FrameStep = 1;
DefaultJob.Parameters.Images.ColorChannel = 1;

% Output (vector) parameters
DefaultJob.Parameters.Vectors.Directory = '~/Desktop/piv_images/vect';
DefaultJob.Parameters.Vectors.BaseName = 'frame_';
DefaultJob.Parameters.Vectors.NumberOfDigits = 6;

% Grid parameters
% These are the same for all processing methods
DefaultJob.Parameters.Processing.Grid.Spacing.X = 32; % Horizontal spacing between grid points (pixels)
DefaultJob.Parameters.Processing.Grid.Spacing.Y = 32; % Vertical spacing between grid points (pixels)
DefaultJob.Parameters.Processing.Grid.Buffer.Y = [32 32]; % Top and bottom grid buffer (pixels)
DefaultJob.Parameters.Processing.Grid.Buffer.X = [32 32]; % Left and right grid buffer (pixels)

% Interrogation Region Parameters
% These are the same for all processing methods
DefaultJob.Parameters.Processing.InterrogationRegion.Height = 128;
DefaultJob.Parameters.Processing.InterrogationRegion.Width = 128;
DefaultJob.Parameters.Processing.InterrogationRegion.SpatialWindowFraction = [0.5 0.5];
DefaultJob.Parameters.Processing.InterrogationRegion.ZeroMeanRegion = 1;

% FMC options
DefaultJob.Parameters.Processing.InterrogationRegion.FMIWindowSize = [2, 2, 0];
DefaultJob.Parameters.Processing.InterrogationRegion.FMIWindowType = 'hann1';

% Correlation parameters
% Correlation type
DefaultJob.Parameters.Processing.Correlation.Method = 'rpc'; 

% Peak fit method ('least-squares', ...)
DefaultJob.Parameters.Processing.Correlation.PeakFitMethod =...
    '3-point';

% FMC parameters
% FMC Image resampling parameters
DefaultJob.Parameters.Processing.Correlation.FMC.NumberOfRings = 64;
DefaultJob.Parameters.Processing.Correlation.FMC.NumberOfWedges = 256;
DefaultJob.Parameters.Processing.Correlation.FMC.MinimumRadius = 2;

% FMC parameters for windowing the FMI log polar images.
DefaultJob.Parameters.Processing.Correlation.FMC.FMIWindowSize = [2, 2, 0];
DefaultJob.Parameters.Processing.Correlation.FMC.FMIWindowType = 'hann1';


% RPC filter diameter for the spatial correlation (pixels)
DefaultJob.Parameters.Processing.Correlation.FMC.SpatialRPCDiameter = 2.8; 

% RPC spectral filter diameter for the FMC correlation (pixels)
DefaultJob.Parameters.Processing.Correlation.FMC.FMCDiameter = 3.3; 

DefaultJob.Parameters.Processing.Correlation.FMC.FMCFilterType = 'relative';

% FFT Dimensions.
% Is this specific to FMC? 
DefaultJob.Parameters.Processing.FFTSize = [128, 128];

% Search multiple peaks?
DefaultJob.Parameters.Processing.MultiPeak = 0;

% Get rid of this?
DefaultJob.Parameters.Processing.FMC.FmcDifferenceMethod = 'forward';

% Discrete window offset parameters
DefaultJob.Parameters.Processing.DWO.DoDiscreteWindowOffset = 0;
DefaultJob.Parameters.Processing.DWO.DwoDifferenceMethod = 'central';

% Image deformation parameters
DefaultJob.Parameters.Processing.Deform.DoImageDeformation = 0;
DefaultJob.Parameters.Processing.Deform.MaxIterations = 1;
DefaultJob.Parameters.Processing.Deform.ConvergenceCriterion = 0.01;

% Smoothing parameters
DefaultJob.Parameters.Processing.Smoothing.DoSmoothing = 0;
DefaultJob.Parameters.Processing.Smoothing.KernelDiameter = 7; 
DefaultJob.Parameters.Processing.Smoothing.KernelGaussianStdDev = 1;

% Universal Outlier Detection Parameters
DefaultJob.Parameters.Processing.Validation.DoValidation = 1;
DefaultJob.Parameters.Processing.Validation.UodStencilRadius = 1;
DefaultJob.Parameters.Processing.Validation.UodThreshold = 3;
DefaultJob.Parameters.Processing.Validation.UodExpectedDifference = [0.1, 0.1];
DefaultJob.Parameters.Processing.Validation.UodMedianThreshold = [3, 2];
DefaultJob.Parameters.Processing.Validation.UodWindowSize = [3, 3; 3, 3];
DefaultJob.Parameters.Processing.Validation.UThresh = [-inf, inf];
DefaultJob.Parameters.Processing.Validation.VThresh = [-inf, inf];

% Default Processing parameters
defaultProcessing = DefaultJob.Parameters.Processing;



% Job 1, Pass 1
SegmentItem = DefaultJob;
SegmentItem.Parameters.Images.CorrelationStep = 1;
SegmentItem.JobOptions.NumberOfPasses = 1;

SegmentItem.JobOptions.StartFromExistingField = 0;
SegmentItem.JobOptions.StartPass = 1;

% Pass 1
SegmentItem.Parameters.Processing(1) = defaultProcessing;

% Get rid of this
SegmentItem.Parameters.Processing(1).DwoDifferenceMethod = 'central';
%

SegmentItem.Parameters.Processing(1).Grid.Spacing.X = 64;
SegmentItem.Parameters.Processing(1).Grid.Spacing.Y = 64;
SegmentItem.Parameters.Processing(1).Grid.Buffer.Y = [0, 0];
SegmentItem.Parameters.Processing(1).Grid.Buffer.X = [0, 0];
SegmentItem.Parameters.Processing(1).InterrogationRegion.Height = 128;
SegmentItem.Parameters.Processing(1).InterrogationRegion.Width = 128;
SegmentItem.Parameters.Processing(1).DWO.DoDiscreteWindowOffset = 1;
SegmentItem.Parameters.Processing(1).Deform.DoImageDeformation = 0;
SegmentItem.Parameters.Processing(1).Smoothing.DoSmoothing = 0;
SegmentItem.Parameters.Processing(1).Correlation.Method = 'rpc';
SegmentItem.Parameters.Processing(1). ...
    InterrogationRegion.SpatialWindowFraction = [0.50 0.50];

SegmentItem.Parameters.Processing(2) = SegmentItem.Parameters.Processing(1);

% Add to job list
JOBLIST(1) = SegmentItem;



end

