%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
addpath('Data_Import_Functions', 'Other_Useful_Functions');

tic_beginning = tic;

%% SET VALUES

% DataFolder_path = 'G:/My Drive/EWS-MultipleWindow-Method/Data';
DataFolder_path = '../Data';

System_name = 'PowerSystem';

time_transient = 20;
overlap_ratio = 80/100;
smallest_step_size = 20000;

% Set significance values
significance_value_tau = 0.05;
gpu_shift_critical_size = 520;

H_val_to_match = -1;


% Define function to call to print parfor progress bar
progress_bar_data_queue = parallel.pool.DataQueue;
afterEach(progress_bar_data_queue, @Progress_Bar_func);


%% IMPORT DATA AND PREPROCESS

% Import data into 'Data' struct
Data = Import_Data(DataFolder_path, System_name);

time = Data.time;
parameter_variation = Data.parameter_variation;
state_timeseries = Data.state_timeseries;

parameter_bifurcation = Data.parameter_bifurcation;
bifurcation_time = Data.bifurcation_time;
rate_of_parameter_variation = Data.rate_of_parameter_variation;

time_transient;
sampling_frequency = Data.sampling_frequency;
delta_t = Data.delta_t;

% Remove Data struct to save space
clearvars Data

% Remove initial transients and plot new time series
selection_transient = time > time_transient;
time = time(selection_transient);
parameter_variation = parameter_variation(selection_transient);
state_timeseries = state_timeseries(selection_transient);

% Plot transient removed state timeseries
figure('Name', 'Timeseries_TransientRemoved');
tiledlayout(1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
nexttile;
plot(time, state_timeseries);
xlabel('Time');
ylabel('State Variable Time Series');
nexttile;
plot(time, parameter_variation);
xlabel('Time');
ylabel('Parameter Time Variation');

fprintf('LIST OF USEFUL DATA\n');
fprintf('--------------------\n');
fprintf('time                   = refer figure\n');
fprintf('parameter_variation    = refer figure\n');
fprintf('state_timeseries       = refer figure\n');
fprintf('parameter_bifurcation  = %f\n', parameter_bifurcation);
fprintf('bifurcation_time       = %f s\n', bifurcation_time);
fprintf('time_transient         = %f s\n', time_transient);
fprintf('sampling_frequency     = %f Hz\n', sampling_frequency);
fprintf('delta_t                = %f s\n', delta_t);
fprintf('\n\n');


%% SET WINDOW DETAILS

overlap_ratio;
max_overlap_ratio = 99 / 100;
min_overlap_ratio = 80 / 100;

largest_window_size = floor(bifurcation_time / delta_t);
largest_step_size = floor(largest_window_size * (1 - overlap_ratio));

smallest_step_size;
smallest_window_size = ceil(smallest_step_size / (1 - max_overlap_ratio));

window_size_increment = smallest_window_size / 10;
window_size_list = smallest_window_size: window_size_increment: largest_window_size;
step_size_list = floor(window_size_list * (1 - overlap_ratio));

total_window_count = length(window_size_list);

fprintf('SET WINDOW DETAILS\n');
fprintf('--------------------\n');
fprintf('largest_window_size   = %10d / %f s\n', largest_window_size, largest_window_size * delta_t);
fprintf('largest_step_size     = %10d / %f s\n', largest_step_size, largest_step_size * delta_t);
fprintf('smallest_window_size  = %10d / %f s\n', smallest_window_size, smallest_window_size * delta_t);
fprintf('smallest_step_size    = %10d / %f s\n', smallest_step_size, smallest_step_size * delta_t);
fprintf('total_window_count    = %10d', total_window_count);
fprintf('\n\n');


%% FIND EWS TIMESERIES FOR EACH WINDOW SIZE

fprintf('EWS TIME SERIES\n');
fprintf('--------------------\n');

% Display progress bar and set progress to 0
Progress_Bar_func('begin', total_window_count);

parfor k = 1: total_window_count

    % Update progress bar each time a new loop starts
    send(progress_bar_data_queue, 'ongoing');
    
    % Set window size and step size to calculate the corresponding EWS timeseries
    window_size = window_size_list(k);
    window_step = step_size_list(k);

    % Generate EWS timeseries for particular window size
    [EWS_details{k}, RateRMS_details{k}] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);
    
end

Progress_Bar_func('ended');

% Plot figure that shows some representative EWS time series evenly spread across window sizes
figure('Name', 'EWS_Representative');
EWS_representative_window_count = 9;
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for k = floor(linspace(1, total_window_count, EWS_representative_window_count) )
    nexttile;
    plot(EWS_details{k}.time_window_ends, EWS_details{k}.AC_timeseries);
    xlim([smallest_window_size * delta_t, time(end)]);
end


%% EXAMINING SIGNIFICANCE OF EWS TIME SERIES TRENDS

fprintf('SIGNIFICANCE VALUES OF EWS TIME SERIES\n');
fprintf('--------------------\n');

Progress_Bar_func('begin', total_window_count);

EWS_timeseries_lengths = zeros(1, total_window_count);

parfor k = 1: total_window_count

    % Update progress bar each time a new loop starts
    send(progress_bar_data_queue, 'ongoing');

    n_ktau = length(EWS_details{k}.time_window_ends);
    EWS_timeseries_lengths(k) = n_ktau;

    % Generate Kendall-tau value and significance level time series. Calculate for increasing amounts of EWS time series data.
    % Preallocate [tau, z, p, H] vectors for increasing amounts EWS time series data
    time_EWS{k} = EWS_details{k}.time_window_ends;
    tau{k} = zeros(1, n_ktau);
    z{k} = zeros(1, n_ktau);
    p{k} = zeros(1, n_ktau);
    H{k} = zeros(1, n_ktau);

    % Setting evaluation parameters
    min_ktau_length = 5;
    significance_value_ac = 0.05;

    for j = min_ktau_length: n_ktau

        % Prepare timeseries to be fed for calculating kendall-tau
        time_ktau = EWS_details{k}.time_window_ends(1: j);
        AC_timeseries_ktau = EWS_details{k}.AC_timeseries(1: j);
        
        % Set significance value level to keep or discard autocorrelation at particular lags amongst the data.
        print_bool = 0;

        % Calculate Kendall-tau and determine whether to reject or retain null hypothesis
        [tau{k}(j), z{k}(j), p{k}(j), H{k}(j)] = Modified_MannKendall_test(time_ktau, AC_timeseries_ktau, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);

    end

end

Progress_Bar_func('ended');

% Plot figure that shows some representative EWS time series evenly spread across window sizes. Also, plot points of maximum significance.
figure('Name', 'EWS_Representative_MaximumSignificance');
EWS_representative_window_count = 9;
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact');

for k = floor(linspace(1, total_window_count, EWS_representative_window_count) )
    nexttile;
    plot(time_EWS{k}, p{k});
    xlim([smallest_window_size * delta_t, time(end)]);
end


%% MAKE PREDICTION MAP

% Normalize window size with largest window size, i.e. the one till bifurcation.
% Normalize time with bifurcation time.

tic_method_1_timer_1 = tic;
% METHOD - 1
% Join time series
my_time = horzcat(time_EWS{:});
my_H = horzcat(H{:});
y_val = repelem((window_size_list / largest_window_size), EWS_timeseries_lengths);

tic_sorting = tic;
% Sort the arrays
[my_time, time_order] = sort(my_time);
my_H = my_H(time_order);
y_val = y_val(time_order);
time_sorting = toc(tic_sorting);

my_H_1 = (my_H == H_val_to_match);
my_H_2 = (my_H == 2);
my_H_0 = ~(my_H_1 | my_H_2);

% Plot the prediction map
figure('Name', 'Prediction_Map_1');

plot(my_time(my_H_1), y_val(my_H_1), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey4, 'MarkerEdgeColor' , PS.Grey5);
plot(my_time(my_H_2), y_val(my_H_2), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
plot(my_time(my_H_0), y_val(my_H_0), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey1, 'MarkerEdgeColor' , PS.Grey2);

xlabel('Normalized Time');
ylabel('Normalized Window Size');

time_method_1_timer_1 = toc(tic_method_1_timer_1);
fprintf('time_method_1_timer_1 = %f\n', time_method_1_timer_1);


tic_method_2_timer_1 = tic;
% METHOD - 2
% Preallocate cell arrays
t_1 = cell(1, total_window_count); t_2 = t_1; t_0 = t_1;
y_val_1 = cell(1, total_window_count); y_val_2 = y_val_1; y_val_0 = y_val_1;

for k = 1: total_window_count

    y_val = window_size_list(k) / largest_window_size;
    
    % Note that below H_1 is for H_val_to_match, which may correspond to H = -1 or 1 as appropriate
    H_1 = (H{k} == H_val_to_match);
    H_2 = (H{k} == 2);
    H_0 = ~(H_1 | H_2);
    
    t_1{k} = time_EWS{k}(H_1) / bifurcation_time;
    t_2{k} = time_EWS{k}(H_2) / bifurcation_time;
    t_0{k} = time_EWS{k}(H_0) / bifurcation_time;

    y_val_1{k} = y_val * ones(1, length(t_1{k}));
    y_val_2{k} = y_val * ones(1, length(t_2{k}));
    y_val_0{k} = y_val * ones(1, length(t_0{k}));
    
end

% Convert cell arrays into vectors to plot
t_1 = horzcat(t_1{:});
t_2 = horzcat(t_2{:});
t_0 = horzcat(t_0{:});

y_val_1 = horzcat(y_val_1{:});
y_val_2 = horzcat(y_val_2{:});
y_val_0 = horzcat(y_val_0{:});

% Plot the prediction map
figure('Name', 'Prediction_Map_2');

plot(t_1, y_val_1, 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey4, 'MarkerEdgeColor' , PS.Grey5);
plot(t_2, y_val_2, 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
plot(t_0, y_val_0, 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey1, 'MarkerEdgeColor' , PS.Grey2);

xlabel('Normalized Time');
ylabel('Normalized Window Size');

time_method_2_timer_1 = toc(tic_method_2_timer_1);
fprintf('time_method_2_timer_1 = %f\n\n', time_method_2_timer_1);


%% PLOT PREDICTION FRACTION FROM SEARCHING ACTUAL DATA

% Prediction fraction at time t = (Number of H values in favor of the trend) / (Total H values) where both are measured at time t




tic_method_1_timer_2 = tic;
% Find and sort the unique time values at which prediction fraction will be calculated
time_prediction_frac = unique(my_time);
% Method - 1
for i = 1: length(time_prediction_frac)

    t = time_prediction_frac(i);

    % Time interval for doing averaging of prediction fraction - here taken to be a multiple of the smallest step size
    % The largest value in diff(time_prediction_frac) = smallest_step_size
    % Do only left side averaging otherwise we will end up using future data to determine current prediction fraction
    n_steps_avg = 10;
    t_interval = n_steps_avg * smallest_step_size * delta_t;
    
    t_val_idx = ( (my_time >= (t - t_interval)) & (my_time <= t) );
    H_taken = my_H(t_val_idx);
    H_total(i) = length(H_taken);
    H_favor(i) = sum(H_taken == H_val_to_match);

end

% Calculate prediction fraction
prediction_fraction_1 = H_favor ./ H_total;

time_method_1_timer_2 = toc(tic_method_1_timer_2);
fprintf('time_method_1_timer_2 = %f\n', time_method_1_timer_2);


tic_method_2_timer_2 = tic;
% Find and sort the unique time values at which prediction fraction will be calculated
time_prediction_frac = unique(horzcat(time_EWS{:}));
% Method - 2

% Preallocate vectors to hold number of H values in favor and the total H values at time t
H_favor = zeros(1, length(time_prediction_frac));
H_total = zeros(1, length(time_prediction_frac));

% Calculate values of the above variables at each instance
for i = 1: length(time_prediction_frac)

    t = time_prediction_frac(i);

    % Time interval for doing averaging of prediction fraction - here taken to be a multiple of the smallest step size
    % The largest value in diff(time_prediction_frac) = smallest_step_size
    % Do only left side averaging otherwise we will end up using future data to determine current prediction fraction
    n_steps_avg = 10;
    t_interval = n_steps_avg * smallest_step_size * delta_t;

    % Check if this time is available for each EWS timeseries
    for k = 1: total_window_count
        for t_val = time_prediction_frac( (time_prediction_frac >= (t - t_interval)) & (time_prediction_frac <= t) )
            % If available then add H values to corresponding vectors
            t_idx = find(time_EWS{k} == t_val);
            if ~isempty(t_idx)
                H_total(i) = H_total(i) + 1;
                if H{k}(t_idx) == H_val_to_match
                    H_favor(i) = H_favor(i) + 1;
                end
            end
        end

    end

end

% Calculate prediction fraction
prediction_fraction = H_favor ./ H_total;

time_method_2_timer_2 = toc(tic_method_2_timer_2);
fprintf('time_method_2_timer_2 = %f\n\n', time_method_2_timer_2);



figure('Name', 'Prediction_Fraction_SearchMethod');
% Plot the prediction fraction vs normalized time
plot(time_prediction_frac / bifurcation_time, prediction_fraction, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Normalized Time');
ylabel('Prediction Fraction');




fprintf('time_sorting = %f\n', time_sorting);


total_time_of_running = toc(tic_beginning);
fprintf('total_time_of_running = %f\n\n', total_time_of_running);



%% SAVE ALL FIGURES

figure_folder_name = sprintf('RPV_%.5f__OR_%.3f__SL_%.6f__SW_%d__TT_%.3f', rate_of_parameter_variation, overlap_ratio, significance_value_tau, smallest_window_size, time_transient);
figure_location = sprintf('Prediction_Maps/%s/%s', System_name, figure_folder_name);

total_figure_count = length(findobj('type','figure'));

% Check if the directory already exists, if it doesn't, create it
if not(isfolder(figure_location))
    mkdir(figure_location)
end

% Save the figures
for i = 1: total_figure_count
    fig_id = figure(i);
    savefig(fig_id, sprintf('%s/%s.fig', figure_location, fig_id.Name));

    % If figure is prediction map save it as bmp
    if strcmp(fig_id.Name, 'Prediction_Map_2')
        pred_map_image_file_location = sprintf('%s/%s.png', figure_location, fig_id.Name);
        exportgraphics(fig_id, pred_map_image_file_location, 'Resolution', 1200);
    end

end


return

%% PLOT PREDICTION FRACTION FROM IMAGE OF PREDICTION MAP

Colors_used.color_yes = [];
Colors_used.color_no = [];

prediction_fraction_from_Image = PredictionFraction_from_Image_func(pred_map_image_file_location, Colors_used);

















