%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
addpath('Data_Import_Functions');


%% SET VALUES

DataFolder_path = 'G:/My Drive/EWS-MultipleWindow-Method/Data';

System_name = 'PowerSystem';

time_transient = 20;


%% IMPORT DATA AND PREPROCESS

% Import data into 'Data' struct
Data = Import_Data(DataFolder_path, System_name);

time = Data.time;
parameter_variation = Data.parameter_variation;
state_timeseries = Data.state_timeseries;

parameter_bifurcation = Data.parameter_bifurcation;
bifurcation_time = Data.bifurcation_time;

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

fprintf('LIST OF USEFUL TIMESERIES\n');
fprintf('--------------------\n');
fprintf('time                   = refer figure\n');
fprintf('parameter_variation    = refer figure\n');
fprintf('state_timeseries       = refer figure\n');
fprintf('parameter_bifurcation  = %f\n', parameter_bifurcation);
fprintf('bifurcation_time       = %f s\n', bifurcation_time);
fprintf('sampling_frequency     = %f Hz\n', sampling_frequency);
fprintf('delta_t                = %f s\n', delta_t);
fprintf('\n\n');











