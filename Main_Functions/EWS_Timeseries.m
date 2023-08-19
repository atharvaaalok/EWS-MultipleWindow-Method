function [EWS_timeseries, time_window_ends] = EWS_Timeseries(time, state_timeseries, window_size, window_step, chosen_EWS)

%% REMOVE TAIL DATA

% Remove data at the end that will be left after the last window
data_total_length = length(time);
tail_data_remove_length = mod((data_total_length - window_size), window_step);

time = time(1: end - tail_data_remove_length);
state_timeseries = state_timeseries(1: end - tail_data_remove_length);


%% CALCULATE EWS MEASURES

window_idx = 1: window_step: (length(time) - window_size + 1);

% Get the time values at which the EWS values will be generated. The EWS will be represented at the ends of the windows
time_window_ends = time(window_idx + window_size - 1);

% Generate the specified EWS time series
switch chosen_EWS
    case 'Skewness'
        EWS_timeseries = rolling_skewness(state_timeseries, window_size, window_step);
    case 'Kurtosis'
        EWS_timeseries = rolling_kurtosis(state_timeseries, window_size, window_step);
    case 'Autocorrelation'
        EWS_timeseries = rolling_autocorr_lag1(state_timeseries, window_size, window_step);
    otherwise
        warning('Unsupported EWS measure.');
        exit();
end


end