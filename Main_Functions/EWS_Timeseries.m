function [EWS_details, RateRMS_details] = EWS_Timeseries(time, state_timeseries, parameter_variation, delta_t, window_size, window_step)

%% CONVERT VECTORS TO GPU ARRAYS

tic

% Convert vector to GPU array for much faster computations

gpu_available = canUseGPU();
if gpu_available == 1
    time = gpuArray(time);
    state_timeseries = gpuArray(state_timeseries);
end

a = toc;


%% REMOVE TAIL DATA

tic

data_total_length = length(time);
tail_data_remove_length = mod((data_total_length - window_size), window_step);

time = time(1: end - tail_data_remove_length);
state_timeseries = state_timeseries(1: end - tail_data_remove_length);
parameter_variation = parameter_variation(1: end - tail_data_remove_length);

% fprintf('REMOVE TAIL DATA\n');
% fprintf('--------------------\n');
% fprintf('data_total_length \t\t\t\t\t= %d\n', data_total_length);
% fprintf('tail_data_remove_length \t\t\t= %d\n', tail_data_remove_length);
% fprintf('tail_data_remove_length in seconds \t= %.3f s\n', tail_data_remove_length * delta_t);
% fprintf('\n\n');

b = toc;


%% SET WINDOW INDICES AND TIME VECTOR OF WINDOW ENDS

tic

window_idx = 1: window_step: (length(time) - window_size + 1);

time_window_ends_idx = window_idx(1:end) + window_size - 1;

time_window_ends = time(time_window_ends_idx);
parameter_window_ends = parameter_variation(time_window_ends_idx);

% fprintf('SET WINDOW INDICES AND TIME VECTOR OF WINDOW ENDS\n');
% fprintf('--------------------\n');
% fprintf('data_total_length - window_idx(end) \t= %d\n', data_total_length - window_idx(end));
% fprintf('time_window_ends(end) \t\t\t\t\t= %.3f s\n', time_window_ends(end));
% fprintf('time(end) \t\t\t\t\t\t\t\t= %.3f s\n', time(end));
% fprintf('\n\n');

c = toc;


%% ALLOCATE SPACE FOR EWS MEASURES AND KENDALL-TAU VALUES AND SET PROPERTIES

tic

% RMS
rms_timeseries = zeros(1, length(time_window_ends));

% Autocorrelation
AC_lag = 1;
AC_timeseries = zeros(1, length(time_window_ends));

% Note: The delta_t input to the function is useful when AC_lag is to be set in seconds

d = toc;


%% CALCULATE EWS MEASURES

tic

for i_window = 1: length(window_idx)
    timeseries_window_data = state_timeseries(window_idx(i_window): window_idx(i_window) + window_size - 1);
    
    % RMS
    rms_timeseries(i_window) = rms(timeseries_window_data);

    % Autocorrelation
%     acf = autocorr(timeseries_window_data, NumLags = AC_lag);
%     AC_timeseries(i_window) = acf(AC_lag + 1);
    AC_timeseries(i_window) = skewness(timeseries_window_data);

end

e = toc;


%% RATE OF CHANGE OF RMS

tic

drms_state = diff(rms_timeseries);
dt = diff(time_window_ends);
rate_rms_state = drms_state ./ dt;

[max_rate_rms_state, max_rate_rms_state_idx] = max(rate_rms_state);
time_vector_rate = time_window_ends(2: end);
time_max_rate = time_vector_rate(max_rate_rms_state_idx);

parameter_vector_rate = parameter_window_ends(2: end);
parameter_max_rate = parameter_vector_rate(max_rate_rms_state_idx);

rate_transition_time = time_max_rate;

% fprintf('RATE OF CHANGE OF RMS\n');
% fprintf('--------------------\n');
% fprintf('max_rate_rms_pressure \t= %.2f Pa/s\n', max_rate_rms_state);
% fprintf('time_max_rate \t\t\t= %.2f s\n', time_max_rate);
% fprintf('parameter_max_rate \t\t= %.2f\n', parameter_max_rate);
% fprintf('\n\n');

f = toc;


%% SET RETURN VALUES

tic

if gpu_available == 1
    % If some vectors had been made GPU array then convert them back to normal arrays
    
    % Make EWS timeseries object to return all timeseries succinctly
    EWS_details.time_window_ends = gather(time_window_ends);
    EWS_details.parameter_window_ends = gather(parameter_window_ends);
    EWS_details.rms_timeseries = gather(rms_timeseries);
    EWS_details.AC_timeseries = gather(AC_timeseries);
    
    % Return rate of rms details as an object for succinctness
    RateRMS_details.rate_rms_state = gather(rate_rms_state);
    RateRMS_details.max_rate_rms_state = gather(max_rate_rms_state);
    RateRMS_details.time_vector_rate = gather(time_vector_rate);
    RateRMS_details.time_max_rate = gather(time_max_rate);
    RateRMS_details.parameter_vector_rate = gather(parameter_vector_rate);
    RateRMS_details.parameter_max_rate = gather(parameter_max_rate);
    RateRMS_details.rate_transition_time = gather(rate_transition_time);
else
    % Make EWS timeseries object to return all timeseries succinctly
    EWS_details.time_window_ends = time_window_ends;
    EWS_details.parameter_window_ends = parameter_window_ends;
    EWS_details.rms_timeseries = rms_timeseries;
    EWS_details.AC_timeseries = AC_timeseries;
    
    % Return rate of rms details as an object for succinctness
    RateRMS_details.rate_rms_state = rate_rms_state;
    RateRMS_details.max_rate_rms_state = max_rate_rms_state;
    RateRMS_details.time_vector_rate = time_vector_rate;
    RateRMS_details.time_max_rate = time_max_rate;
    RateRMS_details.parameter_vector_rate = parameter_vector_rate;
    RateRMS_details.parameter_max_rate = parameter_max_rate;
    RateRMS_details.rate_transition_time = rate_transition_time;
end


g = toc;

TS_gen_total_time = a + b + c + d + e + f + g;

% fprintf("a = %f\n", a);
% fprintf("b = %f\n", b);
% fprintf("c = %f\n", c);
% fprintf("d = %f\n", d);
% fprintf("e = %f\n", e);
% fprintf("f = %f\n", f);
% fprintf("g = %f\n", g);
% fprintf("TS_gen_total_time = %f\n", TS_gen_total_time);


end