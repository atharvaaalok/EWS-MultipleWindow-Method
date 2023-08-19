function AC_timeseries = rolling_autocorr_lag1(x, window_size, window_step)
    
    sum_x = [0, cumsum(x)];
    sum_x2 = [0, cumsum(x.^2)];
    sum_xy = [0, cumsum(x(2: end) .* x(1: end - 1))];
    
    % Get the indices at which to calculate autocorrelation lag 1 values
    window_ends_idx = window_size: window_step: length(x);
    
    % Preallocate array for storing autocorrelation lag 1 time series
    AC_timeseries = zeros(1, length(window_ends_idx));

    % Generate autocorrelation lag 1 time series
    N = window_size - 1;
    N2 = N * N;
    for i = 1: length(window_ends_idx)
        idx = window_ends_idx(i);
        iwin1 = idx - window_size + 1;
        sum_use = sum_x(idx + 1) - sum_x(iwin1);
        sum_x_lag1 = sum_use - x(idx);
        sum_x_fwd1 = sum_use - x(iwin1);
        sum_xy_window = sum_xy(idx) - sum_xy(iwin1);
        sum_x2_window = sum_x2(idx + 1) - sum_x2(iwin1);
        var_dev_window = (sum_x2_window - sum_use.^2 / window_size) / (window_size - 1); % square of std_dev_window
        % Calculate autocorrelation lag 1 for the window
        AC_timeseries(i) = (N * sum_xy_window - sum_x_lag1 * sum_x_fwd1) / (N2 * var_dev_window);
    end

end