function AC_timeseries = rolling_autocorr_lag1(x, window_size, window_step)
    
    N = length(x);
    
    sum_x = zeros(1, N);
    sum_x(1) = x(1);
    sum_xy = zeros(1, N);

    for n = 2: N
        sum_x(n) = sum_x(n - 1) + x(n);
        sum_xy(n) = sum_xy(n - 1) + x(n) * x(n - 1);
    end
    
    % Get the indices at which to calculate autocorrelation lag 1 values
    window_ends_idx = window_size: window_step: length(x);
    
    % Preallocate array for storing autocorrelation lag 1 time series
    AC_timeseries = zeros(1, length(window_ends_idx));

    % Generate autocorrelation lag 1 time series
    for i = 1: length(window_ends_idx)
        idx = window_ends_idx(i);
        if i == 1
            sum_use = sum_x(idx);
        else
            sum_use = sum_x(idx) - sum_x(idx - window_size);
        end
        sum_x_lag1 = sum_use - x(idx);
        sum_x_fwd1 = sum_use - x(idx - window_size + 1);
        sum_xy_window = sum_xy(idx) - sum_xy(idx - window_size + 1);
        std_dev_window = std(x(idx - window_size + 1: idx));
        N = window_size - 1;
        % Calculate autocorrelation lag 1 for the window
        AC_timeseries(i) = (N * sum_xy_window - sum_x_lag1 * sum_x_fwd1) / (N * std_dev_window * N * std_dev_window);
    end

end