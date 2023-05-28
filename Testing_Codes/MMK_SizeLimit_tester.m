%% Modified Mann-Kendall Test GPU SWITCH CRITICAL LENGTH
% The Modified Mann-Kendall test function converts normal arrays into GPU arrays to speed up the computations for calculating median in the Sen slope part of the code
% But when the arrays are too short then the overhead time to transfer the array to GPU memory takes more time than the median calculation for a normal array
% This code finds the minimum length above which using the GPU array is beneficial

clear; clc; close all;
PS = PLOT_STANDARDS();


%% COMPUTE THE TIME FOR 

gpu_use_list = [0, 1];

GPU_time = [];
NoGPU_time = [];

for k = 1: length(gpu_use_list)
    
    % Set whether to use the GPU or not
    gpu_use = gpu_use_list(k);

    n = 2000;
    significance_value_tau = 0.05;
    significance_value_ac = 0.05;
    print_bool = 0;

    time_taken = [];
    
    for i = 5: n
        % i
        X = rand(1, i);
        time = 1: i;
        [tau, z, p, H, time_taken] = MMK(time, X, significance_value_tau, significance_value_ac, gpu_use, print_bool, time_taken);
    end
    
    time_taken = time_taken';
    
    % Save the time taken for each length of vector in corresponding array depending on whether GPU was used or not.
    if gpu_use == 1
        GPU_time = time_taken;
    else
        NoGPU_time = time_taken;
    end

end


%% PLOT THE VARIATIONS OF TIME TAKEN IN CASE OF 1) GPU USAGE 2) NO GPU USAGE

N = min([length(GPU_time), length(NoGPU_time)]);
GPU_time = GPU_time(1: N);
NoGPU_time = NoGPU_time(1: N);

% Plot the variation of time taken for GPU and No GPU use cases
figure(1);
hold on
plot(1: N, GPU_time, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor', PS.Red1);
plot(1: N, NoGPU_time, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Blue2, 'MarkerEdgeColor', PS.Blue2);

% Plot the variation for shorter array lengths to find where the GPU computation speed overtakes normal array computations
limit_length = 1000;
difference_time = GPU_time(1: limit_length) - NoGPU_time(1: limit_length);
difference_time = difference_time';

k = floor((5/100) * limit_length);
[low_remove, low_remove_idx] = mink(difference_time, k);
[large_remove, large_remove_idx] = maxk(difference_time, k);

selection = ones(1, limit_length) == 1;
selection(low_remove_idx) = 0;
selection(large_remove_idx) = 0;

difference_time = difference_time(selection);

figure(2);
plot(1: length(difference_time), difference_time, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Green1, 'MarkerEdgeColor', PS.Green1);
ylim([min(difference_time), max(difference_time)]);









