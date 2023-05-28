clear; clc; close all;
PS = PLOT_STANDARDS();


GPU_time = load('GPU_MMK_time.mat');
GPU_time = GPU_time.var;

NoGPU_time = load('NoGPU_MMK_time.mat');
NoGPU_time = NoGPU_time.var;
% 
% figure(1);
% hold on
% plot(1: length(GPU_time), GPU_time, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Red1, 'MarkerEdgeColor', PS.Red1);
% plot(1: length(NoGPU_time), NoGPU_time, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Blue2, 'MarkerEdgeColor', PS.Blue2);

difference_time = GPU_time - NoGPU_time(1: end - 2);
vec = difference_time(10: 1000);
t = 1: length(vec);
figure(2);
plot(t, vec, 'LineStyle', 'none', 'Marker', '.', 'MarkerFaceColor', PS.Green1, 'MarkerEdgeColor', PS.Green1);
ylim([min(vec), max(vec)]);