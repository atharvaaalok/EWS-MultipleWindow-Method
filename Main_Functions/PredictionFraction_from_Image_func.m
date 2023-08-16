function prediction_fraction = PredictionFraction_from_Image_func(image_file_location, Colors_used)
    
    % Import Image
    Color_mat = imread(image_file_location);
    Color_mat = double(Color_mat);

    % Scaling array for image array dimensionality reduction
    scaling_array = [256^2, 256, 1];

    % Reduce dimensions of image array and color arrays
    color_no_1D = Colors_used.color_no * scaling_array;
    color_yes_1D = Colors_used.color_yes * scaling_array;
    Color_mat_1D = reshape(Color_mat, [], 3) * scaling_array;

    % Find the pixel locations of color match and reshape back to 2D
    pixels_favor = reshape(ismember(Color_mat_1D, color_yes_1D), size(Color_mat_1D, 1), []);
    pixels_not_favor = reshape(ismember(Color_mat_1D, color_no_1D), size(Color_mat_1D, 1), []);
    
    N_favor = sum(pixels_favor);
    N_not_favor = sum(pixels_not_favor);

    prediction_fraction = N_favor ./ (N_favor + N_not_favor);

end