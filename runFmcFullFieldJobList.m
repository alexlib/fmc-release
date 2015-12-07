function runFmcFullFieldJobList(JOBLIST)

% Add compiled code path
addpath ba_interpolation;

% Determine the number of jobs in the job list.
nJobs = length(JOBLIST);

% Determine the best fft algoritm to use. 
fftw('planner', 'patient');
    
% Run each job in the job list
for n = 1 : nJobs
    
    % Start a timer.
    job_tic = tic;
    
    % Extract the jobfile from the joblist.
    JobFile = JOBLIST(n);
    
    % Grid parameters
    gridSpacingX = JobFile.Parameters.Processing(1).Grid.Spacing.X;
    gridSpacingY = JobFile.Parameters.Processing(1).Grid.Spacing.Y;
   
    % Interrogation Region Parameters
    region_height = ...
        JobFile.Parameters.Processing(1).InterrogationRegion.Height;
    region_width = ...
        JobFile.Parameters.Processing(1).InterrogationRegion.Width; 
    
    % Correlation method (just read here for naming purposes)
    correlation_method = JobFile.Parameters. ...
        Processing(1).Correlation.Method;
    
    % Determine whether to skip image pairs whose results exist.
    skipExisting = JobFile.JobOptions.SkipExisting;
    
    % File path information for the input images
    image_directory = JobFile.Parameters.Images.Directory;
    image_base_name = JobFile.Parameters.Images.BaseName;
    image_extension = JobFile.Parameters.Images.Extension;
    
    % Number of digits in the output file names.
    image_number_of_digits = JobFile.Parameters.Images.NumberOfDigits;
    start_image = JobFile.Parameters.Images.Start;
    end_image = JobFile.Parameters.Images.End;
    frame_step = JobFile.Parameters.Images.FrameStep;
    correlation_step = JobFile.Parameters.Images.CorrelationStep;
    
    % File path information for the output vector fields
    output_directory = JobFile.Parameters.Vectors.Directory;
    output_base_name = JobFile.Parameters.Vectors.BaseName;
    
    % Number of digits in the output file names.
    output_number_of_digits = JobFile.Parameters.Vectors.NumberOfDigits;
    
    % Default output base name in case the field
    % was empty in the jobfile.
    % This block specifies the output base name
    % according to the processing parameters
    % of the first PIV pass.
    if isempty(output_base_name)
        [image_base_name correlation_method ...
            '_grid' num2str(gridSpacingX) 'x' ...
            num2str(gridSpacingY) '_region' ...
            num2str(region_height) 'x' ...
            num2str(region_width) '_'];
    end
        
    % Default output directory in case the field
    % was empty in the jobfile
    % This block specifies the output directory 
    % to be next to the image directory
    % and names it according to the processing parameters
    % of the first PIV pass.
    if isempty(output_directory)
        output_directory = fullfile(...
            image_directory, '..', 'vect', ...
            correlation_method, ...
            [num2str(region_height) 'x' num2str(region_width)],...
            ['cstep_' num2str(correlation_step, '%02.0f')]);        
    end
    
    % Determine the path to the output directory and create it if it doesn't exist.
    if ~exist(output_directory, 'dir')
        mkdir(output_directory);
    end
    
    
    % Image sets
    % (delete)
%     startSet = JobFile.Parameters.Sets.Start;
%     endSet   = JobFile.Parameters.Sets.End;
%     setVector = startSet : endSet;
%     nSets = length(setVector);
  
    % Loop over the sets.
    for s = 1 : nSets;
    
    % Strings specifying the number format for the images and output files
    image_number_format  = ['%0' num2str(image_number_of_digits)  '.0f'];
    output_number_format = ['%0' num2str(output_number_of_digits) '.0f'];
    
% Run the job

    % Build a list of image numbers
    firstImageNumbers = start_image : frame_step : end_image;
    secondImageNumbers = firstImageNumbers + correlation_step;
    
    % Determine the number of images
    nPairs = length(firstImageNumbers);
    
    % Build a list of image file paths
    for k = 1 : nPairs
        
        firstImageFilePaths(k, :) = fullfile(image_directory, [image_base_name num2str(firstImageNumbers(k), image_number_format) image_extension]);
        secondImageFilePaths(k, :) = fullfile(image_directory, [image_base_name num2str(secondImageNumbers(k), image_number_format) image_extension]);
        outputFilePath(k, :) = fullfile(output_directory, [output_base_name num2str(firstImageNumbers(k), image_number_format) '_' num2str(secondImageNumbers(k), image_number_format) '.mat']);
        
        FilePaths(k).FirstImagePath = firstImageFilePaths(k, :);
        FilePaths(k).SecondImagePath = secondImageFilePaths(k, :);
        FilePaths(k).OutputFilePath = outputFilePath(k, :);
        
    end

    % Do the correlations
    for k = 1 : nPairs
        
        if ~(skipExisting && exist(FilePaths(k).OutputFilePath, 'file'))
            disp(['Correlating pair ' num2str(k) ' of ' num2str(nPairs)]);
            pair_tic = tic;
            fmcFullField(FilePaths(k), JobFile);
            disp(['Saved vector field to ' outputFilePath(k, :)]);
            disp(['Image Pair Time: ' num2str(toc(pair_tic)) ' sec'])
        end
    end
   
    
    end % End looping over sets   

    job_toc = toc(job_tic);
    fprintf('Total job time: %d seconds\n', job_toc);
    
end


end


