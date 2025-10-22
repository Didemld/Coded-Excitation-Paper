clc
clear all
close all
%%
initialize_k
addpath('D:\OneDrive - Delft University of Technology\1TUD\Athesis\All_codes\CRLB\Data_128\128x5')
%% Source = burst
tone_burst_freq = 15*1e6;          % [Hz]
tone_burst_cycles = 3;
sampling_freq = 1/kgrid.dt;     % [Hz] 
steering_angle = 0;            % [deg]
element_spacing = dx;           % [m]
base_pulse = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);
numelements = num_elements;
%% encoding matrix C
pulse_length = length(base_pulse);
code_length = 5;
%%
load A_tr_128x5_120&128_200&65_down4.mat
A_real = real(A_tr);
[V,D] = eig(A_real);  % calculate eigenvectors and eigenvalues of A
[lambda, idx] = max(diag(D));  % find largest eigenvalue and its index
%% 
a = flip(diag(D));
k = 10*log10(a);
figure
plot(real(k(1:360)),'LineWidth',3);
ylabel('10*log10(Eigenvalue)')
title('Eigenvalues of A_{tr}')
xlim tight