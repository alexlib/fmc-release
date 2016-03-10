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
Job.Parameters.Images.Directory = '/Users/matthewgiarra/Desktop/test_images/raw';
Job.Parameters.Images.BaseName = 'Img_';
Job.Parameters.Images.Extension = '.tif';
Job.Parameters.Images.NumberOfDigits = 3;

% Output file parameters parameters
Job.Parameters.Vectors.Directory = '/Users/matthewgiarra/Desktop/test_images/vect';
Job.Parameters.Vectors.BaseName = 'fmc_proc2_';
Job.Parameters.Vectors.NumberOfDigits = 3;

% Start and end images
Job.Parameters.Images.Start = 15;
Job.Parameters.Images.End = 15;
Job.Parameters.Images.FrameStep = 1;
Job.Parameters.Images.CorrelationStep = 1;

% Masking
Job.Parameters.Mask.DoMasking = false;
%Job.Parameters.Mask.Path = '~/Desktop/mask.tif';

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
Job.Parameters.Processing(1).Correlation.Method = 'fmc';

% Grid and region parameters.
% Region height and width refer to the un-windowed
% interrogation region, i.e., double what Prana 
% calls the "window resolution."
Job.Parameters.Processing(1).InterrogationRegion.Height = 96;
Job.Parameters.Processing(1).InterrogationRegion.Width = 96;

% Grid spacing
Job.Parameters.Processing(1).Grid.Spacing.X = 24;
Job.Parameters.Processing(1).Grid.Spacing.Y = 24;

% Post processing
Job.Parameters.Processing(1).Smoothing.DoSmoothing = false;
Job.Parameters.Processing(1).Validation.DoValidation = false;

% Iterative scheme
Job.Parameters.Processing.Iterative.Method = 'none';
Job.Parameters.Processing(1).Iterative.MaxIterations = 1;

% % Make the second pass!
% Copy the parameters from the first pass to a new pass.
Job.Parameters.Processing(2) = Job.Parameters.Processing(1);

% Modify some parameters for the second pass.
Job.Parameters.Processing(2).Grid.Spacing.X = 24;
Job.Parameters.Processing(2).Grid.Spacing.Y = 24;
Job.Parameters.Processing(2).InterrogationRegion.Height = 96;
Job.Parameters.Processing(2).InterrogationRegion.Width  = 96;

% Add the segment to job list
JOBLIST(1) = Job;

% % Add a second job to the job list
% % Change the processing methods to FMC
% Job.Parameters.Processing(1).Correlation.Method = 'fmc';
% Job.Parameters.Processing(2).Correlation.Method = 'fmc';
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

