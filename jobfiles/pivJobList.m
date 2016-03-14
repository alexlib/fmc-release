function JOBLIST = pivJobList()

% This line specifies whether you want to run the 
% list of jobs after creating it (i.e., during the call to
% pivJobList). If true, calling pivJobList with no arguments
% will create the list of jobs specified in this file and then 
% run them in order. If false, it will just create the 
% list of jobs, and return them as a structure called JOBLIST.
run_job_list = true;

% Load the defaults from a separate file.
% This file contains all of the parameters
% necessary to run a job, but many of them 
% can stand to be hidden from the user
% in most cases. Feel free to tinker with it 
% and swee what happens.
Job = DefaultJob;

% Default processing parameters. 
% This is useful if you want to do a bunch
% of jobs that all use the same (or similar) processing
% and you don't want to modify everything every time.
default_processing = Job.Parameters.Processing(1);

% Input file parameters
Job.Parameters.Images.Directory = '~/Documents/School/VT/Research/Aether/schlieren/analysis/data/jhu_buyoant_turbulence/raw';
Job.Parameters.Images.BaseName = 'jhu_buoy_turb_';
Job.Parameters.Images.Extension = '.tif';
Job.Parameters.Images.NumberOfDigits = 2;

% Output file parameters parameters
Job.Parameters.Vectors.Directory = '~/Documents/School/VT/Research/Aether/schlieren/analysis/data/jhu_buyoant_turbulence/vect';
Job.Parameters.Vectors.BaseName = 'frame_';
Job.Parameters.Vectors.NumberOfDigits = 2;

% Start and end images
Job.Parameters.Images.Start = 1;
Job.Parameters.Images.End = 1;
Job.Parameters.Images.FrameStep = 1;
Job.Parameters.Images.CorrelationStep = 1;

% Masking
Job.Parameters.Mask.DoMasking = false;
Job.Parameters.Mask.Path = '~/Desktop/mask.tif';

% Number of passes for this job.
% This will override the number
% of passes contained in the jobfile; it's done
% this way so that you can quickly truncate
% the number of passes to perform without deleting
% any text from this file.
Job.JobOptions.NumberOfPasses = 1;

% Processing parameters for the first pass.
Job.Parameters.Processing(1) = default_processing;

% Correlation type to use
Job.Parameters.Processing(1).Correlation.Method = 'scc';

% Grid and region parameters.
% Region height and width refer to the un-windowed
% interrogation region, i.e., double what Prana 
% calls the "window resolution."
Job.Parameters.Processing(1).InterrogationRegion.Height = 128;
Job.Parameters.Processing(1).InterrogationRegion.Width = 128;

% Grid spacing
Job.Parameters.Processing(1).Grid.Spacing.X = 16;
Job.Parameters.Processing(1).Grid.Spacing.Y = 16;
Job.Parameters.Processing(1).Grid.Buffer.Y = [8, 8];
Job.Parameters.Processing(1).Grid.Buffer.X = [8, 8];

% Post processing
Job.Parameters.Processing(1).Smoothing.DoSmoothing = true;
Job.Parameters.Processing(1).Validation.DoValidation = true;

% Iterative scheme
Job.Parameters.Processing(1).Iterative.Method = 'none';
Job.Parameters.Processing(1).Iterative.MaxIterations = 3;

% Attempt to converge the iterative method?
Job.Parameters.Processing(1).Iterative.Converge = true;

% % Make the second pass!
% Copy the parameters from the first pass to a new pass.
Job.Parameters.Processing(2) = Job.Parameters.Processing(1);

% Iterative method options
Job.Parameters.Processing(2).Iterative.MaxIterations = 6;
Job.Parameters.Processing(2).Iterative.Converge = true;

% Modify some parameters for the second pass.
Job.Parameters.Processing(2).Grid.Spacing.X = 8;
Job.Parameters.Processing(2).Grid.Spacing.Y = 8;
Job.Parameters.Processing(2).InterrogationRegion.Height = 64;
Job.Parameters.Processing(2).InterrogationRegion.Width  = 64;

% Add the segment to job list
JOBLIST(1) = Job;

% Correlation type to use
Job.Parameters.Processing(1).Correlation.Method = 'rpc';
JOBLIST(end + 1) = Job;

% % Add a second job to the job list
% % Change the processing methods to FMC
% Job.Parameters.Processing(1).Correlation.Method = 'rpc';
% Job.Parameters.Processing(2).Correlation.Method = 'rpc';
% 
% 
% 
% % Append the job to the job list.
% JOBLIST(end + 1) = Job;

% This line runs the jobfile.
% The option run_job_list is set
% at the very top of this file.
if run_job_list
    runPivFullFieldJobList(JOBLIST);
end

end

