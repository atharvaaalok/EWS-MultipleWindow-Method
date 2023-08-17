%% INITIAL SETUP

clear; clc; close all;
addpath(genpath('Main_Functions'));
% Import useful data and method parameters
Useful_Data_for_Method()

tic_beginning = tic;


%% SET VALUES

tic_set_values = tic;

% DataFolder_path = 'G:/My Drive/EWS-MultipleWindow-Method/Data';
DataFolder_path = '../Data';

System_name = 'PowerSystem';

time_transient = 20;
overlap_ratio = 99/100;
smallest_step_size = 20000;

% Set significance values
significance_value_tau = 0.05;
gpu_shift_critical_size = 520;

H_val_to_match = 1;

% Define function to call to print parfor progress bar
progress_bar_data_queue = parallel.pool.DataQueue;
afterEach(progress_bar_data_queue, @Progress_Bar_func);

fprintf('Time Taken: Set Values = %f\n\n', toc(tic_set_values));


%% IMPORT DATA AND PREPROCESS

tic_import_data = tic;

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

transient_select_num = find(time >= time_transient, 1);
time = time(transient_select_num: end);
parameter_variation = parameter_variation(transient_select_num: end);
state_timeseries = state_timeseries(transient_select_num: end);

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

fprintf('Time Taken: Import Data = %f\n\n', toc(tic_import_data));


%% SET WINDOW DETAILS

tic_window_details = tic;

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

fprintf('Time Taken: Window Details = %f\n\n', toc(tic_window_details));


%% FIND EWS TIMESERIES FOR EACH WINDOW SIZE

tic_EWS_timeseries = tic;

fprintf('EWS TIME SERIES\n');
fprintf('--------------------\n');

% Display progress bar and set progress to 0
Progress_Bar_func('begin', total_window_count);

for k = 1: total_window_count

    % Set window size and step size to calculate the corresponding EWS timeseries
    window_size = window_size_list(k);
    window_step = step_size_list(k);

    % Generate EWS timeseries for particular window size
    [EWS_details{k}, RateRMS_details{k}] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);

    % Update progress bar each time a new loop starts
    send(progress_bar_data_queue, 'ongoing');
    
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

fprintf('Time Taken: Window Details = %f\n\n', toc(tic_EWS_timeseries));

return


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

return
%% MAKE PREDICTION MAP: NORMALIZED WINDOW SIZE AND NORMALIZED TIME

tic_pred_map = tic;

% Normalize window size with largest window size, i.e. the one till bifurcation.
% Normalize time with bifurcation time.

% Join arrays together to extract favorable and unfavorable regions using vector operations
t_EWS_all = horzcat(time_EWS{:});
H_EWS_all = horzcat(H{:});
y_EWS_all = repelem((window_size_list / largest_window_size), EWS_timeseries_lengths);

% Sort time array and 
[t_EWS_all, t_EWS_order_idx] = sort(t_EWS_all);
H_EWS_all = H_EWS_all(t_EWS_order_idx);
y_EWS_all = y_EWS_all(t_EWS_order_idx);

H_1 = (H_EWS_all == H_val_to_match);
H_2 = (H_EWS_all == 2);
H_0 = ~(H_1 | H_2);

% Plot the prediction map
figure('Name', 'Prediction_Map');
hold on
plot(t_EWS_all(H_2), y_EWS_all(H_2), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 1.5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
plot(t_EWS_all(H_0), y_EWS_all(H_0), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 1.5, 'MarkerFaceColor', PS.Grey2, 'MarkerEdgeColor' , PS.Grey2);
plot(t_EWS_all(H_1), y_EWS_all(H_1), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 1.5, 'MarkerFaceColor', PS.Grey5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Normalized Time');
ylabel('Normalized Window Size');

fprintf('Time Taken: Prediction Map = %f\n\n', toc(tic_pred_map));


%% EVALUATE PREDICTION FRACTION FROM ACTUAL DATA: THE SEARCH METHOD

tic_pred_frac = tic;

% Prediction fraction at time t = (Number of H values in favor of the trend) / (Total H values) where both are measured at time t

% Get all the unique time values at which prediction fraction will be calculated
time_prediction_frac = unique(t_EWS_all);

% Preallocate arrays for storing total windows predicting at time t and total in favor out of those
H_total = zeros(1, length(time_prediction_frac));
H_favor = zeros(1, length(time_prediction_frac));

for i = 1: length(time_prediction_frac)

    t = time_prediction_frac(i);

    % Time interval for doing averaging of prediction fraction - here taken to be a multiple of the smallest step size
    % The largest value in diff(time_prediction_frac) is = smallest_step_size
    % Do only left side averaging otherwise we will end up using future data to determine current prediction fraction
    n_steps_avg = 10;
    t_interval = n_steps_avg * smallest_step_size * delta_t;

    t_interval_idx = ( (t_EWS_all >= (t - t_interval)) & (t_EWS_all <= t) );
    H_interval = H_EWS_all(t_interval_idx);
    H_total(i) = length(H_interval);
    H_favor(i) = sum(H_interval == H_val_to_match);

end

% Calculate prediction fraction
prediction_fraction = H_favor ./ H_total;

figure('Name', 'Prediction_Fraction_SearchMethod');
% Plot the prediction fraction vs normalized time
plot(time_prediction_frac / bifurcation_time, prediction_fraction, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Normalized Time');
ylabel('Prediction Fraction');

fprintf('Time Taken: Prediction Fraction = %f\n\n', toc(tic_pred_frac));


%% SAVE ALL FIGURES

tic_save_figures = tic;

figure_folder_name = sprintf('RPV_%.5f__OR_%.3f__SL_%.6f__SW_%d__TT_%.3f', rate_of_parameter_variation, overlap_ratio, significance_value_tau, smallest_window_size, time_transient);
figure_location = sprintf('Prediction_Maps/%s/%s', System_name, figure_folder_name);

% Delete the directory if it already exists. This is to prevent any possible confusion with files left over from a previous run which are not overwritten.
if isfolder(figure_location)
    rmdir(figure_location, 's');
end
% Create the directory
mkdir(figure_location);

total_figure_count = length(findobj('type','figure'));

% Save the figures
for i = 1: total_figure_count
    fig_id = figure(i);
    savefig(fig_id, sprintf('%s/%s.fig', figure_location, fig_id.Name));

    % If figure is prediction map save it as bmp
    if strcmp(fig_id.Name, 'Prediction_Map')
        pred_map_image_file_location = sprintf('%s/%s.png', figure_location, fig_id.Name);
        exportgraphics(fig_id, pred_map_image_file_location, 'Resolution', 1200);
    end

end

fprintf('Time Taken: Save Figures = %f\n\n', toc(tic_save_figures));


%% EVALUATE PREDICTION FRACTION FROM IMAGE OF PREDICTION MAP: THE IMAGE METHOD

tic_pred_frac_image = tic;

% Set color values used for image generation
Colors_used.color_yes = PS.Grey5;
Colors_used.color_no = PS.Grey2;

prediction_fraction_from_Image = PredictionFraction_from_Image_func(pred_map_image_file_location, Colors_used);
time_pred_frac_from_Image = linspace(t_EWS_all(1), t_EWS_all(end), length(prediction_fraction_from_Image));

figure('Name', 'Prediction_Fraction_ImageMethod');
% Plot the prediction fraction vs normalized time
plot(time_pred_frac_from_Image / bifurcation_time, prediction_fraction_from_Image, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Normalized Time');
ylabel('Prediction Fraction');

fprintf('Time Taken: Prediction Fraction from Image = %f\n\n', toc(tic_pred_frac_image));



total_time_of_running = toc(tic_beginning);
fprintf('total_time_of_running = %f\n\n', total_time_of_running);










