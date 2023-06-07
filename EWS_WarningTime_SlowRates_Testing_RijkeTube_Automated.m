%% AUTOMATING PREDICTION MAP GENERATION

% This file automates the generation of prediction maps for different values of the following parameters:
% 1) Time transient
% 2) Overlap ratio
% 3) Smallest step size
% 4) Significance value for tau


%% INITIAL SETUP

clear; clc; close all;


%% DEFINE ALL THE CONDITIONS FOR WHICH PREDICTION MAPS HAVE TO BE CREATED

% The parameter value list is in the following order:
% 1) Time Transient, 2) Overlap ratio, 3) Smallest step size, 4) Significance value for tau
Parameter_list = [
% OR - 80
0, 0.80, 1500, 0.05;
0, 0.80, 1500, 0.01;
0, 0.80, 1500, 0.005;
0, 0.80, 1500, 0.001;

0, 0.80, 1000, 0.05;

0, 0.80, 500, 0.05;
0, 0.80, 500, 0.01;
0, 0.80, 500, 0.005;
0, 0.80, 500, 0.001;

% OR - 90
0, 0.90, 1500, 0.05;

0, 0.90, 1000, 0.05;

0, 0.90, 500, 0.05;

% OR - 95
0, 0.95, 1500, 0.05;

0, 0.95, 1000, 0.05;

0, 0.95, 500, 0.05;

% OR - 99
0, 0.99, 1500, 0.05;
0, 0.99, 1500, 0.01;
0, 0.99, 1500, 0.005;
0, 0.99, 1500, 0.001;
0, 0.99, 1500, 0.0001;
0, 0.99, 1500, 0.00001;

0, 0.99, 1000, 0.05;
0, 0.99, 1000, 0.01;
0, 0.99, 1000, 0.005;
0, 0.99, 1000, 0.001;

0, 0.99, 500, 0.05;
0, 0.99, 500, 0.01;
0, 0.99, 500, 0.005;
0, 0.99, 500, 0.001;
0, 0.99, 500, 0.0001;
0, 0.99, 500, 0.00001;
];

% Create table with appropriate variable headings
Parameter_table = array2table(Parameter_list, "VariableNames", ["time_transient", "overlap_ratio", "smallest_step_size", "significance_value_tau"]);

for i = 1: height(Parameter_table)

    % Set the parameter values to be used for generating the prediction map
    time_transient = Parameter_table.time_transient(i);
    overlap_ratio = Parameter_table.overlap_ratio(i);
    smallest_step_size = Parameter_table.smallest_step_size(i);
    significance_value_tau = Parameter_table.significance_value_tau(i);

    % Generate Prediction maps for a given set of parameter values
    EWS_WarningTime_SlowRates_Testing_PowerSystem_func(time_transient, overlap_ratio, smallest_step_size, significance_value_tau);

    % Close all generated plots
    close all;

end





%% CODE FOR THE PREDICITION MAP GENERATION FUNCTION

function EWS_WarningTime_SlowRates_Testing_PowerSystem_func(time_transient, overlap_ratio, smallest_step_size, significance_value_tau)
%% INITIAL SETUP

% clear; clc; close all;
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

% Set transient time limit
time_transient;

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
fprintf('time_transient \t\t\t\t\t= %.2f\n', time_transient);
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


%% REMOVE INTIAL TRANSIENTS AND PLOT NEW TIMESERIES

time_transient = time_transient;
selection_transient = time > time_transient;
time = time(selection_transient);
state_timeseries = state_timeseries(selection_transient);
parameter_variation = parameter_variation(selection_transient);

% Plot transient removed state timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on
plot(time, state_timeseries);


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
fprintf('largest_window_size \t\t\t\t= %d\n', largest_window_size);
fprintf('largest_step_size \t\t\t\t\t= %d\n', largest_step_size);
fprintf('largest window size in seconds \t\t= %fs\n', largest_window_size * delta_t);
fprintf('largest step size in seconds \t\t= %fs\n', largest_step_size * delta_t);
fprintf('smallest_window_size \t\t\t\t= %d\n', smallest_window_size);
fprintf('smallest_step_size \t\t\t\t\t= %d\n', smallest_step_size);
fprintf('smallest window size in seconds \t= %fs\n', smallest_window_size * delta_t);
fprintf('smallest step size in seconds \t\t= %fs\n', smallest_step_size * delta_t);
fprintf('\n\n');


%% FIND EWS TIMESERIES FOR EACH WINDOW-SIZE

loop_start_tic = tic;

for k = 1: length(window_size_list)

    k
    
    % Set window size and step size to calculate the corresponding EWS timeseries
    window_size = window_size_list(k);
    window_step = step_size_list(k);
    window_length_in_seconds = window_size * delta_t;
    window_step_in_seconds = window_step * delta_t;
    
    EWS_tic = tic;

    % Generate EWS timeseries for particular window size
    [EWS_details, RateRMS_details] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);
    
    t1 = toc(EWS_tic);

    % Plot RMS timeseries on top of state timeseries
    figure(10);
    hold on
    plot(time, state_timeseries);
    plot(EWS_details.time_window_ends, EWS_details.rms_timeseries );
    % Plot AC timeseries
    figure(11);
    plot(EWS_details.time_window_ends, EWS_details.AC_timeseries );
    
    tau_tic = tic;

    % Generate Kendall-tau timeseries. Calculate for increasing amount of EWS timeseries data
    % Prepare [tau, z, p, H] vectors for increasing EWS timeseries data
    n_ktau = length(EWS_details.time_window_ends);
    time_EWS{k} = EWS_details.time_window_ends;
    tau{k} = zeros(1, n_ktau);
    z{k} = zeros(1, n_ktau);
    p{k} = zeros(1, n_ktau);
    H{k} = zeros(1, n_ktau);
    for j = 5: n_ktau
        [k, j];
        % Prepare timeseries to be fed for calculating kendall-tau
        time_ktau = EWS_details.time_window_ends(1: j);
        AC_timeseries_ktau = EWS_details.AC_timeseries(1: j);
        
        % Set significance values
        significance_value_tau;
        significance_value_ac = 0.05;
        gpu_shift_critical_size = 520;

        if j == n_ktau
            print_bool = 1;
        else
            print_bool = 0;
        end

        % Calculate Kendall-tau and determine whether to reject or retain null hypothesis
        [tau{k}(j), z{k}(j), p{k}(j), H{k}(j)] = Modified_MannKendall_test(time_ktau, AC_timeseries_ktau, significance_value_tau, significance_value_ac, gpu_shift_critical_size, print_bool);

    end

    t2 = toc(tau_tic);
    window_total_time = t1 + t2;

    fprintf("t1 = %f\n", t1);
    fprintf("t2 = %f\n", t2);
    fprintf("window_total_time = %f\n", window_total_time);

end

loop_time = toc(loop_start_tic);
fprintf("loop_time = %f\n", loop_time);


%% PLOT THE PREDICTION MAP WITH Y-AXIS AS: WINDOW COUNT

figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

y_val = 0;

for k = 1: length(window_size_list)
    
    y_val = y_val + 1;
    
    H_1 = (H{k} == -1);
    H_2 = (H{k} == 2);
    H_0 = ~(H_1 | H_2);
    
    t_1 = time_EWS{k}(H_1);
    t_2 = time_EWS{k}(H_2);
    t_0 = time_EWS{k}(H_0);
    
    plot(t_1, y_val * ones(1, length(t_1)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey4, 'MarkerEdgeColor' , PS.Grey5);
    plot(t_2, y_val * ones(1, length(t_2)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
    plot(t_0, y_val * ones(1, length(t_0)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey1, 'MarkerEdgeColor' , PS.Grey2);
    
end

xlabel('Time');
ylabel('Window Count');

% Save the figure
prediction_map_windowcount_figure_name = sprintf("Prediction_Maps/RijkeTube/Window_Count/Prediction_Map_WindowCount_TT%d_SW%.2f_OR%.2f_SL%.5f_RPV%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, rate_of_parameter_variation);
saveas(gcf, prediction_map_windowcount_figure_name);


%% PLOT THE PREDICTION MAP WITH Y-AXIS AS: WINDOW SIZE

figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

y_val = min(window_size_list);

for k = 1: length(window_size_list)
    
    y_val = window_size_list(k);
    
    H_1 = (H{k} == -1);
    H_2 = (H{k} == 2);
    H_0 = ~(H_1 | H_2);
    
    t_1 = time_EWS{k}(H_1);
    t_2 = time_EWS{k}(H_2);
    t_0 = time_EWS{k}(H_0);
    
    plot(t_1, y_val * ones(1, length(t_1)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey4, 'MarkerEdgeColor' , PS.Grey5);
    plot(t_2, y_val * ones(1, length(t_2)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
    plot(t_0, y_val * ones(1, length(t_0)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey1, 'MarkerEdgeColor' , PS.Grey2);
    
end

xlabel('Time');
ylabel('Window Size');

% Save the figure
prediction_map_windowsize_figure_name = sprintf("Prediction_Maps/RijkeTube/Window_Size/Prediction_Map_WindowSize_TT%d_SW%.2f_OR%.2f_SL%.5f_RPV%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, rate_of_parameter_variation);
saveas(gcf, prediction_map_windowsize_figure_name);


%% PLOT THE PREDICTION MAP WITH Y-AXIS AS NORMALIZED WINDOW SIZE

figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

% Normalize window size with largest window size, i.e. the one till bifurcation and normalize time with bifurcation time.
largest_window_size;
bifurcation_time;

y_val = min(window_size_list) / largest_window_size;

for k = 1: length(window_size_list)
    
    y_val = window_size_list(k) / largest_window_size;
    
    H_1 = (H{k} == -1);
    H_2 = (H{k} == 2);
    H_0 = ~(H_1 | H_2);
    
    t_1 = time_EWS{k}(H_1) / bifurcation_time;
    t_2 = time_EWS{k}(H_2) / bifurcation_time;
    t_0 = time_EWS{k}(H_0) / bifurcation_time;
    
    plot(t_1, y_val * ones(1, length(t_1)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey4, 'MarkerEdgeColor' , PS.Grey5);
    plot(t_2, y_val * ones(1, length(t_2)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'x', 'MarkerSize', 5, 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor' , PS.Red2);
    plot(t_0, y_val * ones(1, length(t_0)), 'LineStyle', 'none', 'LineWidth', 1, 'Marker', 'o', 'MarkerSize', 5, 'MarkerFaceColor', PS.Grey1, 'MarkerEdgeColor' , PS.Grey2);
    
end

xlabel('Normalized Time');
ylabel('Normalized Window Size');

% Save the figure
prediction_map_windowsizenormalized_figure_name = sprintf("Prediction_Maps/RijkeTube/Window_Size_Normalized/Prediction_Map_WindowSizeNormalized_TT%d_SW%.2f_OR%.2f_SL%.5f_RPV%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, rate_of_parameter_variation);
saveas(gcf, prediction_map_windowsizenormalized_figure_name);


%% PLOT THE PREDICTION FRACTION EVOLUTION WITH TIME: MAKE WITH TIME AND NORMALIZED TIME

% Prediction fraction at time t = (Number of H values in favor of the trend) / (Total H values) both measured at time t

% Find and sort the unique time values at which prediction fraction will be calculated
time_prediction_frac = [];
for k = 1: length(window_size_list)
    time_prediction_frac = [time_prediction_frac, time_EWS{k}];
end

time_prediction_frac = unique(time_prediction_frac);

% Create vectors to hold number of H values in favor and the total H values at time t
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
    for k = 1: length(window_size_list)
        for t_val = time_prediction_frac( (time_prediction_frac >= t - t_interval) & (time_prediction_frac <= t) )
            % If available then add H values to corresponding vectors
            t_idx = find(time_EWS{k} == t_val);
            if ~isempty(t_idx)
                H_total(i) = H_total(i) + 1;
                if H{k}(t_idx) == -1
                    H_favor(i) = H_favor(i) + 1;
                end
            end
        end

    end

end

% Calculate prediction fraction
prediction_fraction = H_favor ./ H_total;

% Plot the prediction fraction vs time
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

plot(time_prediction_frac, prediction_fraction, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Time');
ylabel('Prediction Fraction');

% Save the figure
prediction_frac_figure_name = sprintf("Prediction_Maps/RijkeTube/Prediction_Fraction/Prediction_Fraction_TT%d_SW%.2f_OR%.2f_SL%.5f_RPV%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, rate_of_parameter_variation);
saveas(gcf, prediction_frac_figure_name);

% Plot the prediction fraction vs normalized time
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

plot(time_prediction_frac / bifurcation_time, prediction_fraction, 'LineStyle', 'none', 'Marker', '.', 'MarkerSize', 5, 'MarkerEdgeColor' , PS.Grey5);

xlabel('Normalized Time');
ylabel('Prediction Fraction');

% Save the figure
prediction_frac_normalized_figure_name = sprintf("Prediction_Maps/RijkeTube/Prediction_Fraction_Normalized/Prediction_Fraction_Normalized_TT%d_SW%.2f_OR%.2f_SL%.5f_RPV%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, rate_of_parameter_variation);
saveas(gcf, prediction_frac_normalized_figure_name);




end