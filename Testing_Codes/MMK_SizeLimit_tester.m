clear; clc; close all;

global var
var = [];

n = 5000;

significance_value_tau = 0.05;
significance_value_ac = 0.05;

print_bool = 0;

gpu_use = 1;

for i = 5: n
    i
    X = rand(1, i);
    time = 1: i;
    [tau, z, p, H] = MMK(time, X, significance_value_tau, significance_value_ac, gpu_use, print_bool);
end

var = var';

if gpu_use == 1
    save("GPU_MMK_time.mat", "var");
else
    save("NoGPU_MMK_time.mat", "var");
end