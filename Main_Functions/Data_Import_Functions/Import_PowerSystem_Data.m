function Data = Import_PowerSystem_Data(DataFolder_path)
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
    
    filename = sprintf('NoiseOmega5_mu%.4f_delta%.2f_omega%.2f_E%.2f_Pm%.4f_t%.2f_deltaT%.5f_ConstantTimeStep.mat', mu, delta0, omega0, E0, Pm0, t2, delta_t);
    file_path = sprintf('%s/PowerSystem/Noise5/%s', DataFolder_path, filename);
    load(file_path);
    
    tSol;
    xSol = YSol(:, 1);
    ySol = YSol(:, 2);
    omegaSol = YSol(:, 3);
    ESol = YSol(:, 4);
    PmSol = YSol(:, 5);
    
    % Plot omega timeseries
    figure('Name', 'Timeseries_Original');
    hold on
    plot(tSol, omegaSol);
    xlabel('Time');
    ylabel('$$\\omega$$', 'Interpreter', 'Latex');

    % Parameter bifurcation value
    parameter_bifurcation = 0.6495;
    
    fprintf('IMPORT DATA\n');
    fprintf('--------------------\n');
    fprintf('t_0                          = %f s\n', t1);
    fprintf('t_f                          = %f s\n', t2);
    fprintf('sampling_rate                = %d Hz\n', sampling_rate);
    fprintf('rate_of_parameter_variation  = %f\n', mu);
    fprintf('parameter_0                  = %f\n', PmSol(1));
    fprintf('parameter_f                  = %f\n', PmSol(end));
    fprintf('parameter_bifurcation        = %f\n', parameter_bifurcation);
    fprintf('\n\n');


    %% CONVERT TO GENERAL VARIABLE NAMES AND ASSIGN TO RETURN DATA STRUCT
    
    Data.time = tSol;
    Data.parameter_variation = PmSol;
    Data.state_timeseries = omegaSol;
    
    Data.parameter_bifurcation = parameter_bifurcation;
    Data.rate_of_parameter_variation = mu;
    Data.bifurcation_time = (parameter_bifurcation - Data.parameter_variation(1)) / Data.rate_of_parameter_variation;
    
    Data.sampling_frequency = sampling_rate;
    Data.delta_t = delta_t;



end