function kr_timeseries = rolling_kurtosis(x, window_size, window_step)
    
    % First generate arrays for incremental Moments
    N = length(x);
    M1_list = zeros(1, N);
    M1_list(1) = x(1);
    M2_list = zeros(1, N);
    M3_list = zeros(1, N);
    M4_list = zeros(1, N);

    for n = 2: N
        delta = x(n) - M1_list(n - 1);
        M1_list(n) = M1_list(n - 1) + delta / n;
        M4_list(n) = M4_list(n - 1) + (delta * delta * delta * delta) * (n - 1) * (n * n - 3 * n + 3) / (n * n * n) + 6 * delta * delta * M2_list(n - 1) / (n * n) - 4 * delta * M3_list(n - 1) / n;
        M3_list(n) = M3_list(n - 1) + delta * delta * delta * (n - 1) * (n - 2) / (n * n) - 3 * delta * M2_list(n - 1) / n;
        M2_list(n) = M2_list(n - 1) + delta * delta * (n - 1) / n;
    end
    
    % Get the indices at which to calculate skewness values.
    window_ends_idx = window_size: window_step: length(x);
    
    % Preallocate array for storing skewness time series
    kr_timeseries = zeros(1, length(window_ends_idx));

    % Generate skewness time series
    for i = 1: length(window_ends_idx)
        idx = window_ends_idx(i);
        M2 = M2_list(idx) - M2_list(idx - window_size + 1);
        M4 = M4_list(idx) - M4_list(idx - window_size + 1);
        % Calculate skewness for the window
        kr_timeseries(i) = (M4 * window_size) / (M2 * M2);
    end

end