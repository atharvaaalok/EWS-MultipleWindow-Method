function skewness = skewness_onepass(x)
    
    N = length(x);
    M1 = 0;
    M2 = 0;
    M3 = 0;

    for n = 1: N
        delta = x(n) - M1;
        M1 = M1 + delta / n;
        M3 = M3 + delta * delta * delta * (n - 1) * (n - 2) / (n * n) - 3 * delta * M2 / n;
        M2 = M2 + delta * delta * (n - 1) / n;
    end
    
    std_dev = sqrt(M2 / N);

    % Calculate skewness
    skewness = (M3 / N) / (std_dev * std_dev * std_dev);
    
end