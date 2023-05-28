%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
figure_counter = 0;


%% IMPORT DATA

% Load pressure data from file
data_folder_location = '.';
pressure_file_location_in_data_folder = '../Data/Slow_Rates/Rijke_Tube';
pressure_file_number = 4;
pressure_file_location = sprintf('%s/%s/%d.txt', data_folder_location, pressure_file_location_in_data_folder, pressure_file_number);

Timeseries_data = load(pressure_file_location);
time = Timeseries_data(:, 1);
Pressure_timeseries = Timeseries_data(:, 2);

% Convert to row vectors
time = reshape(time, 1, length(time));
Pressure_timeseries = reshape(Pressure_timeseries, 1, length(Pressure_timeseries));

% Pre-process timeseries data
time = time - time(1);
pressure_conversion_factor = (1000 / 0.2175);
Pressure_timeseries = Pressure_timeseries * pressure_conversion_factor;

% Load voltage, current and power value
control_parameter_in_pressure_file_folder = '.';
% The 0 in the following list is due to 18th file being absent
rate_of_parameter_variation_list = [40, 2, 4, 5, 10, 20, 30, 60, 80, 120, 240, 480, 800, 1200, 2400, 4800, 24000, 0, 3];
rate_of_parameter_variation = rate_of_parameter_variation_list(pressure_file_number);
control_parameter_file_location = sprintf('%s/%s/%s/%d_parameter_%.2fmV_s.txt', data_folder_location, pressure_file_location_in_data_folder, control_parameter_in_pressure_file_folder, pressure_file_number, rate_of_parameter_variation);

% Pre-process parameter variation data
Parameter_variation_data = load(control_parameter_file_location);
parameter_voltage_factor = 1.6;
parameter_current_factor = 80;
Voltage_raw_data = Parameter_variation_data(:, 1) * parameter_voltage_factor;
Current_raw_data = Parameter_variation_data(:, 2) * parameter_current_factor;

% Convert to row vectors
Voltage_raw_data = reshape(Voltage_raw_data, 1, length(Voltage_raw_data));
Current_raw_data = reshape(Current_raw_data, 1, length(Current_raw_data));
Power_raw_data = Voltage_raw_data .* Current_raw_data;

% Sampling frequency
pressure_sampling_frequency_list = [10, 4, 4, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 10, 4, 4] * 1000;
pressure_sampling_frequency = pressure_sampling_frequency_list(pressure_file_number);
parameter_sampling_frequency_list = [2, 1, 1, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 0, 1] * 1000;
parameter_sampling_frequency = parameter_sampling_frequency_list(pressure_file_number);

% Fit line to parameter variation data
raw_data_time = linspace(time(1), time(end), length(Power_raw_data));
[xData, yData] = prepareCurveData(raw_data_time, Power_raw_data);
ft = fittype('poly4');
[fitresult, gof] = fit(xData, yData, ft, 'Normalize', 'on');
Power_data = fitresult(time);

% Fit line to voltage variation data
[xData, yData] = prepareCurveData(raw_data_time, Voltage_raw_data);
ft = fittype('poly1');
[fitresult, gof] = fit(xData, yData, ft);
Voltage_data = fitresult(time);

% Plot pressure timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
plot(time, Pressure_timeseries);
xlabel('Time (s)');
ylabel('Pressure (Pa)');

% Plot pressure with Power
figure_counter = figure_counter + 1;
figure(figure_counter);
plot(Power_data, Pressure_timeseries);
xlabel('Power (J)');
ylabel('Pressure (Pa)');


fprintf('IMPORT DATA\n');
fprintf('--------------------\n');
fprintf('pressure_file_number \t\t\t= %d\n', pressure_file_number);
fprintf('t_0 \t\t\t\t\t\t\t= %.2f s\n', time(1));
fprintf('t_f \t\t\t\t\t\t\t= %.2f s\n', time(end));
fprintf('pressure_sampling_frequency \t= %d Hz\n', pressure_sampling_frequency);
fprintf('parameter_sampling_frequency \t= %d Hz\n', parameter_sampling_frequency);
fprintf('rate_of_parameter_variation \t= %.2f mV/s\n', ((Voltage_data(end) - Voltage_data(1)) / (time(end) - time(1))) * 1000);
fprintf('parameter_0 \t\t\t\t\t= %.2fJ\n', Power_data(1));
fprintf('parameter_f \t\t\t\t\t= %.2fJ\n', Power_data(end));
fprintf('\n\n');


%% CONVERT TO GENERAL VARIABLE NAMES

time = time;
state_timeseries = Pressure_timeseries;
parameter_variation = Power_data;
sampling_frequency = pressure_sampling_frequency;
delta_t = 1 / sampling_frequency;
% Bifurcation point
parameter_bifurcation = 590;
Voltage_bifurcation = Voltage_data(Power_data < parameter_bifurcation);
parameter_bifurcation_voltage = Voltage_bifurcation(end);
rate_of_parameter_variation = (rate_of_parameter_variation * parameter_voltage_factor) / 1000;
bifurcation_time = (parameter_bifurcation_voltage - Voltage_data(1)) / rate_of_parameter_variation;


%% LIST OF USEFUL TIMESERIES AND VARIABLES

fprintf('LIST OF USEFUL TIMESERIES\n');
fprintf('--------------------\n');
fprintf('time \t\t\t\t\t\t\t= time\n');
fprintf('state_timeseries \t\t\t\t= Pressure_timeseries\n');
fprintf('parameter_variation \t\t\t= Power_data\n');
fprintf('sampling_frequency \t\t\t\t= pressure_sampling_frequency\n');
fprintf('delta_t \t\t\t\t\t\t= 1 / sampling_frequency\n');
fprintf('parameter_bifurcation \t\t\t= 590J\n');
fprintf('rate_of_parameter_variation \t= rate_of_parameter_variation / 1000\n');
fprintf('\n\n');


%% SET WINDOW DETAILS

overlap_ratio = 99 / 100;

largest_window_size = floor(bifurcation_time / delta_t);
largest_step_size = floor(largest_window_size * (1 - overlap_ratio));

smallest_step_size = 2;
smallest_window_size = ceil(smallest_step_size / (1 - overlap_ratio));

step_size_increment = 1;
step_size_list = smallest_step_size: step_size_increment: largest_step_size;
window_size_list = floor(step_size_list / (1 - overlap_ratio));

fprintf('SET WINDOW DETAILS\n');
fprintf('--------------------\n');
fprintf('largest_window_size \t\t\t\t= %d\n', largest_window_size);
fprintf('largest_step_size \t\t\t\t\t= %d\n', largest_step_size);
fprintf('largest window size in seconds \t\t= %fs\n', largest_window_size * delta_t);
fprintf('largest step size in seconds \t\t= %fs\n', largest_step_size * delta_t);
fprintf('smallest_window_size \t\t\t\t= %d\n', smallest_window_size);
fprintf('smallest_step_size \t\t\t\t\t= %d\n', smallest_step_size);
fprintf('smallest window size in seconds \t= %fs\n', smallest_window_size * delta_t);
fprintf('smallest step size in seconds \t\t= %fs\n', smallest_step_size * delta_t);
fprintf('\n\n');




























