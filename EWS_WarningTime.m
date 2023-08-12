%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
addpath('Data_Import_Functions');


%% SET VALUES

DataFolder_path = 'G:/My Drive/EWS-MultipleWindow-Method/Data';

System_name = 'PowerSystem';

time_transient = 20;
overlap_ratio = 90/100;
smallest_step_size = 40000;


%% IMPORT DATA AND PREPROCESS

% Import data into 'Data' struct
Data = Import_Data(DataFolder_path, System_name);

time = Data.time;
parameter_variation = Data.parameter_variation;
state_timeseries = Data.state_timeseries;

parameter_bifurcation = Data.parameter_bifurcation;
bifurcation_time = Data.bifurcation_time;

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
figure();
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

fprintf('SET WINDOW DETAILS\n');
fprintf('--------------------\n');
fprintf('largest_window_size   = %10d / %f s\n', largest_window_size, largest_window_size * delta_t);
fprintf('largest_step_size     = %10d / %f s\n', largest_step_size, largest_step_size * delta_t);
fprintf('smallest_window_size  = %10d / %f s\n', smallest_window_size, smallest_window_size * delta_t);
fprintf('smallest_step_size    = %10d / %f s\n', smallest_step_size, smallest_step_size * delta_t);
fprintf('\n\n');


%% FIND EWS TIMESERIES FOR EACH WINDOW SIZE

fprintf('EWS TIME SERIES\n');
fprintf('--------------------\n');
fprintf('Current window count/Total window count = ');

figure();
EWS_representative_window_count = 9;
tiledlayout('flow', 'TileSpacing', 'compact', 'Padding', 'compact')

parfor k = 1: length(window_size_list)
    
    % Print the window count for which the loop runs
    printed_length = fprintf('%d/%d', k, length(window_size_list));
    
    % Set window size and step size to calculate the corresponding EWS timeseries
    window_size = window_size_list(k);
    window_step = step_size_list(k);

    % Generate EWS timeseries for particular window size
    [EWS_details{k}, RateRMS_details{k}] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);
    
    % Plot representative EWS time series evenly spread across window sizes
    if ismember(k, floor(linspace(1, length(window_size_list), EWS_representative_window_count) ) )
        nexttile;
        plot(EWS_details{k}.time_window_ends, EWS_details{k}.AC_timeseries);
        xlim([smallest_window_size * delta_t, time(end)]);
    end

    fprintf(repmat('\b', 1, printed_length));

end

fprintf('%d/%d\n', length(window_size_list), length(window_size_list));
fprintf('\n\n');





















