function [x, y, u, v, spectral_planes] = correlate_image_pairs(ImagePaths, PassParameters)

	% Read the file paths
	image_file_paths_01 = ImagePaths.FirstImages;
	image_file_paths_02 = ImagePaths.SecondImages;
	
	% Read the size of one of the images
	% Hopefully the first image exists
	if exist(image_file_paths_01{1})
		[image_height, image_width, num_channels] = size(imread(image_file_paths_01{1}));
	else
		fprintf('Error: failed to load image:\n%s\n', image_file_paths_01{1});
		fprintf('Aborting job.\n')
		return;
	end
	
	% Number of pairs, i.e., the ensemble length
	ensemble_length = length(image_file_paths_01);
	
	% Grid the image
	[x, y] = gridImage([image_height, image_width], ...
            [gridSpacingY gridSpacingX], ...
            gridBufferY, gridBufferX);




end