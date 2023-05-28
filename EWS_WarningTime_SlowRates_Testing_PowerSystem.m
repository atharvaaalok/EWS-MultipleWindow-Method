%% INITIAL SETUP

clear; clc; close all;
PS = PLOT_STANDARDS();
figure_counter = 0;


%% IMPORT DATA

delta0 = 1;
x0 = cos(delta0);
y0 = sin(delta0);
omega0 = 0.9;
E0 = 0.8;
Pm0 = .5;

% Time Range details
% nsteps = 1000000;
sampling_rate = 5001;
delta_t = 1 / (sampling_rate - 1);     % the actual formula should be 1 / (sampling_rate - 1), but I use this as an approximation as integer multiple (5000) makes 1 second.
t1 = 0;

Y0 = [x0; y0; omega0; E0; Pm0];

% mu_list = 0.001:0.0005:0.008;
% t2_list = 300 * ones(1, length(mu_list));

mu = [0.0001];
t2 = 2000;

filename = sprintf('Data/Slow_Rates/Power_System/Rate_0.0001.mat');
load(filename);

tSol;
xSol = YSol(:, 1);
ySol = YSol(:, 2);
omegaSol = YSol(:, 3);
ESol = YSol(:, 4);
PmSol = YSol(:, 5);

% Plot omega timeseries
figure_counter = figure_counter + 1;
figure(figure_counter);
hold on
plot(tSol, omegaSol);

xlabel('Time');
ylabel('$$\\omega$$', 'Interpreter', 'Latex');

% Set transient time limit
time_transient = 20;

fprintf('IMPORT DATA\n');
fprintf('--------------------\n');
fprintf('t_0 \t\t\t\t\t\t\t= %.2f s\n', t1);
fprintf('t_f \t\t\t\t\t\t\t= %.2f s\n', t2);
fprintf('sampling_rate \t\t\t\t\t= %d Hz\n', sampling_rate);
fprintf('rate_of_parameter_variation \t= %.2f mV/s\n', mu);
fprintf('parameter_0 \t\t\t\t\t= %.2f\n', PmSol(1));
fprintf('parameter_f \t\t\t\t\t= %.2f\n', PmSol(end));
fprintf('time_transient \t\t\t\t\t= %.2f\n', time_transient);
fprintf('\n\n');



%% CONVERT TO GENERAL VARIABLE NAMES

time = tSol;
state_timeseries = omegaSol;
parameter_variation = PmSol;
sampling_frequency = sampling_rate;
delta_t = delta_t;
parameter_bifurcation = 0.6495;
rate_of_parameter_variation = mu;
bifurcation_time = (parameter_bifurcation - parameter_variation(1)) / rate_of_parameter_variation;


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
fprintf('time \t\t\t\t\t\t\t= tSol\n');
fprintf('state_timeseries \t\t\t\t= omegaSol\n');
fprintf('parameter_variation \t\t\t= PmSol\n');
fprintf('sampling_frequency \t\t\t\t= sampling_rate\n');
fprintf('delta_t \t\t\t\t\t\t= delta_t\n');
fprintf('parameter_bifurcation \t\t\t= 0.6495\n');
fprintf('rate_of_parameter_variation \t= mu\n');
fprintf('\n\n');


%% SET WINDOW DETAILS

overlap_ratio = 99 / 100;
max_overlap_ratio = 99 / 100;
min_overlap_ratio = 80 / 100;

largest_window_size = floor(bifurcation_time / delta_t);
largest_step_size = floor(largest_window_size * (1 - overlap_ratio));

smallest_step_size = 3000;
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
    
    tic

    % Generate EWS timeseries for particular window size
    [EWS_details, RateRMS_details] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step);
    
    t1 = toc;

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
        significance_value_tau = 0.05;
        significance_value_ac = 0.05;
        gpu_shift_critical_size = 550;

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


%% PLOT THE VARIATION OF P-VALUES AND H

figure_counter = figure_counter + 1;
figure(figure_counter);
hold on

y_val = length(window_size_list);

for k = 1: length(window_size_list)
    
    y_val = y_val - 1;
    
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

xlabel('time');
ylabel('window count');

% Save the figure
prediction_map_figure_name = sprintf("Prediction_Maps/Prediction_Map_TT%d_SW%.2f_OR%.2f_SL%.5f_mu%.4f.fig", time_transient, smallest_window_size, overlap_ratio, significance_value_tau, mu);
saveas(gcf, prediction_map_figure_name);

