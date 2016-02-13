function runPivFullFieldJobList(JOBLIST)
% This function determines which files (images) to correlate
% and the paths to the corresponding output files to be created, 
% and then calls the correlation 

% Add compiled code path
addpath ba_interpolation;

% Determine the number of jobs in the job list.
number_of_jobs = length(JOBLIST);

% Determine the best fft algoritm to use. 
fftw('planner', 'patient');
    
% Run each job in the job list
for n = 1 : number_of_jobs
    
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
    
    % Number of the first image
    start_image = JobFile.Parameters.Images.Start;
    
    % Number of the last image
    end_image = JobFile.Parameters.Images.End;
    
    % Number of images between subsequent frames (frame_step = 1
    % skips no frames).
    frame_step = JobFile.Parameters.Images.FrameStep;
    
    % Correlation step
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
        output_base_name = [image_base_name correlation_method ...
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
      
    % Strings specifying the number format for the images and output files
    image_number_format  = ['%0' num2str(image_number_of_digits)  '.0f'];
    output_number_format = ['%0' num2str(output_number_of_digits) '.0f'];
    
    % Build a list of image numbers
    firstImageNumbers = start_image : frame_step : end_image;
    secondImageNumbers = firstImageNumbers + correlation_step;
    
    % Determine the number of pairs that will be correlated
    number_of_pairs = length(firstImageNumbers);
   
    % Loop over the image pairs.
    for k = 1 : number_of_pairs
        
        % File path to the first image
        FilePaths.FirstImagePath = fullfile(...
            image_directory, [image_base_name...
            num2str(firstImageNumbers(k), image_number_format) ...
            image_extension]);
        
        % File path to the second image
        FilePaths.SecondImagePath = fullfile(...
            image_directory, [image_base_name  ...
            num2str(secondImageNumbers(k), image_number_format) ...
            image_extension]);
        
        % Output file path
        FilePaths.OutputFilePath = fullfile(...
            output_directory, ...
            [output_base_name ...
            num2str(firstImageNumbers(k), output_number_format) ...
            '_' num2str(secondImageNumbers(k), ...
            output_number_format) '_' correlation_method '.mat']);
        
        % This if-statement skips frames for which output data exist
        % if the "skip existing frames" option is enabled.
        if ~(skipExisting && exist(FilePaths.OutputFilePath, 'file'))
            
            % Display a message
            disp(['Correlating pair ' num2str(k) ' of ' num2str(number_of_pairs)]);
            
            % Start a timer for the correlation pair
            pair_tic = tic;
            
            % Perform the correlations on the image pair.
            % This is This is the call to the main code!
            % The try-catch loop assumes that
            % an easy way to break the processing
            % is to specify running the compiled codes
            % but not compiling the codes before hand,
            % so it comiples them if there's an error.
            % % % % % % % % % % % % % % % % % % %
%             try
                pivFullField(FilePaths, JobFile);
%             catch
%                 compile_all
%                 pivFullField(FilePaths, JobFile);
%             end
                
             % % % % % % % % % % % % % % % % % % % 
            
            % Inform the user that the job as completed
            % and print the path to the saved file.
            fprintf('\nSaved vector field to %s\n', ...
                FilePaths.OutputFilePath);
            
            % Inform the user of the total time elapsed for the image pair.
            fprintf('Image Pair Time: %02.f seconds.\n', toc(pair_tic));
        end
        
    end % end of "for k = 1 : number_of_pairs"

    % Calculate and display the elapsed time for the current job.
    job_toc = toc(job_tic);
    fprintf('Total job time: %0.2f seconds\n', job_toc);
    
end % end of "for n = 1 : number_of_jobs"


end % End of function


