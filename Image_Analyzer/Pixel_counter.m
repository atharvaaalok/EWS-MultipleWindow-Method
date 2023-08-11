clear; clc; close all;
PS = PLOT_STANDARDS();



Color_mat = imread("Images/Pred_Map.bmp");

Red_mat = Color_mat(:, :, 1);
Green_mat = Color_mat(:, :, 2);
Blue_mat = Color_mat(:, :, 3);

my_vals = unique(Blue_mat);

my_mat = Red_mat;
my_mat(:, :, 2) = Green_mat;
my_mat(:, :, 3) = Blue_mat;

imwrite(my_mat, "lol.bmp");

flat_mat = zeros(height(my_mat) * width(my_mat), 3);
counter = 0;
for i = 1: height(my_mat)
    for j = 1: width(my_mat)
        counter = counter + 1;
        r = my_mat(i, j, 1);
        g = my_mat(i, j, 2);
        b = my_mat(i, j, 3);
        flat_mat(counter, :) = [r, g, b];
    end
end


[unique_rows, unique_row_idx, similar_idx] = unique(flat_mat, 'rows');

% flat_mat = flat_mat(unique_row_idx, :, :);

count_unique_color = accumarray(similar_idx, 1);


% Get the index of largest 3 values
[new_mat, sort_idx] = sort(count_unique_color);

count_unique_color(sort_idx);
idx_vec = unique_row_idx(sort_idx);

idx_1 = idx_vec(end);
idx_2 = idx_vec(end - 1);
idx_3 = idx_vec(end - 2);

color_1 = flat_mat(idx_1, :);
color_2 = flat_mat(idx_2, :);
color_3 = flat_mat(idx_3, :);

% Add blue in place of all other colors
color_4 = [133, 193, 233];


new_img = zeros(height(my_mat), width(my_mat), 3);

counter = 0;
for i = 1: height(new_img)
    for j = 1: width(new_img)
        counter = counter + 1;
        
        if sum(flat_mat(counter, :) == color_1) == 3
            new_img(i, j, :) = color_1;
        elseif sum(flat_mat(counter, :) == color_2) == 3
%             disp("color 2");
            new_img(i, j, :) = color_2;
        elseif sum(flat_mat(counter, :) == color_3) == 3
%             disp("color 3");
            new_img(i, j, :) = color_3;
        else
            new_img(i, j, :) = color_4;
        end

    end
end



new_img = uint8(new_img);


imwrite(new_img, "New_Pred_Map.bmp");


prediction_fraction = ones(1, width(new_img));

for i = 1: width(new_img)
    grey_counter = 0;
    total_counter = 0;
    for j = 1: height(new_img)
        new_img_color_val = new_img(j, i, :);
        new_img_color_val = new_img_color_val(:);
        new_img_color_val = new_img_color_val';
        if sum(new_img_color_val == color_1) ~= 3
            total_counter = total_counter + 1;
            if sum(new_img_color_val == color_2) == 3
%                 disp("here 2");
                grey_counter = grey_counter + 1;
            end
        end
    end

    prediction_fraction(i) = grey_counter / total_counter;

end


figure(1)
plot(1: length(prediction_fraction), prediction_fraction, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor', PS.Grey5);

saveas(gca, "Pred_factor.png");
