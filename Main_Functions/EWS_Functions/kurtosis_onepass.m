function kurtosis = kurtosis_onepass(x)
    
    N = length(x);
    M1 = 0;
    M2 = 0;
    M3 = 0;
    M4 = 0;

    for n = 1: N
        delta = x(n) - M1;
        M1 = M1 + delta / n;
        M4 = M4 + (delta * delta * delta * delta) * (n - 1) * (n * n - 3 * n + 3) / (n * n * n) + 6 * delta * delta * M2 / (n * n) - 4 * delta * M3 / n;
        M3 = M3 + delta * delta * delta * (n - 1) * (n - 2) / (n * n) - 3 * delta * M2 / n;
        M2 = M2 + delta * delta * (n - 1) / n;
    end

    % Calculate kurtosis
    kurtosis = (M4 * N) / (M2 * M2);
    
end