function autocorr_lag1 = autocorr_lag1_onepass(x)
    
    N = length(x);

    sum_x = x(1);
    sum_xy = 0;

    for n = 2: N
        sum_x = sum_x + x(n);
        sum_xy = sum_xy + x(n) * x(n - 1);
    end
    
    std_dev_x = std(x);
    
    sum_x_lag1 = sum_x - x(N);
    sum_x_fwd1 = sum_x - x(1);

    N = N - 1;

    autocorr_lag1 = (N * sum_xy - sum_x_lag1 * sum_x_fwd1) / (N * std_dev_x * N * std_dev_x);

end