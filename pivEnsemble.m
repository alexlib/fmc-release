function pivEnsemble(JOBFILE)
	
	% Choose whether to run compiled codes.
	run_compiled = JOBFILE.JobOptions.RunCompiled;
	
	% Number of PIV passes
	% Set this to be the lesser of the
	% specified number of passes and the 
	% number of passes for which processing
	% parameters were specified.
	number_of_passes = min(JOBFILE.JobOptions.NumberOfPasses,  ...
	 length(JOBFILE.Parameters.Processing));
	
	 % Extract the processing parameters for 
	 % the number of passes that will be run.
	 JOBFILE.Parameters.Processing = ...
	 	JOBFILE.Parameters.Processing(1 : number_of_passes);

	 % Read the masking parameters
	 do_mask = JOBFILE.Parameters.Mask.DoMasking;
	 
	 % Load the mask if masking is requested
	 if do_mask
		 
		 % Read mask path
		 mask_path = JOBFILE.Parameters.Mask.Path;
		 
		 % Check existence of mask path
		 if exist(mask_path, 'dir');		 
			 % Load the mask
			 mask = double(imread(mask_path));			 
		 else
			 % Inform the user if the mask file wasn't found.
			 fprintf('Error: failed to load the specified mask file:\n%s\n', mask_path);
			 fprintf('Defaulting to no mask.\n');
			 mask = [];
		 end	 
	 end
	
	 % Overwrite the previous pass jobfile and filename info.
	 FilePaths = FILEPATHS;
	 JobFile = JOBFILE;

	 % Initialize the iteration counter 
	 % (this is for DWO and deform convergence)
	 iteration = 1;
	 
	 % Loop over the passes.
	 for pass = 1 : number_of_passes
		 
		 % Read the pass paramteres 
		 PassParameters = JobFile.Parameters.Processing(pass);
		 
		 % Add the mask to the pass parametres
		 PassParameters.Mask = mask;
		 
		 % Image paths
		 ImagePaths = FilePaths.Images;
		 
		 % Loop over the ensemble
		 for pair_num = 1 : ensemble_length
			 
			 % Do the correlations. Note that this
			 % line will work for pair-wise or ensemble
			 % correlations -- the function correlate_image_pair
			 % loops over the lengths of the lists of images
			 % in the ImagePaths structure
			 [x, y, u, v, spectral_planes] = correlate_image_pairs(ImagePaths, PassParameters);
			 
		 end
		 
		 % Increment the pass counter
		 [x, y, u, v, correlation_planes] = piv_image_pair()
		 
	 end
	
end







