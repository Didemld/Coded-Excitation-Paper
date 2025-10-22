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
c_real = V(:, idx);  % extract corresponding complex eigenvector
c_normalized = c_real / norm(c_real);  % normalize to satisfy constraint
C = reshape(c_normalized,numelements,code_length);
Code = C*100;

% C = zeros(numelements,code_length);
% C(:,1) =1;
% C = C./norm(C(:));
% Code =C*100;

% load c_random.mat
% C = temp_array;
% C = C./norm(C(:));
% Code =C*100;
excitation_length = length(base_pulse);
excitation_signals = zeros(excitation_length*code_length,numelements);
for m = 1:numelements
    for k = 1:code_length
       excitation_signals((k-1)*excitation_length+1:k*excitation_length,m) = base_pulse*Code(m,k);
    end
end

%%
filtered_signals=zeros(length(kgrid.t_array),numelements);
for i=1:numelements
    excitation_signal_zeropad = zeros(length(kgrid.t_array),1);
    excitation_signal_zeropad(1:size(excitation_signals,1))=excitation_signals(:,i);
    filtered_signals(:,i) = filterTimeSeries(kgrid, medium, excitation_signal_zeropad);
end
%% Source = burst
% source.p = excitation_signals';
source.p = (filtered_signals)';

%% transform TX to fourier domain

Nt=3500;
TW_zero_padding = zeros(Nt,numelements);
tw_fft = zeros(Nt,1);
TW_fft = zeros(Nt,numelements);

for i = 1:numelements
    TW_zero_padding(:,i) = source.p(i,1:Nt);
    tw_fft =fft(TW_zero_padding(:,i));
    TW_fft(:,i) = [tw_fft(Nt/2+2:end) ;tw_fft(1:Nt/2+1)];
end

fs = 1/kgrid.dt; %sampling frequency
f0 = tone_burst_freq;  %transducer frequency
f_vector = fs*(-(Nt/2-1):(Nt/2))/Nt;

f_min = f0 - f0 *.25;
f_max = f0 + f0 *.25;

f_index = find (f_vector >= f_min & f_vector <= f_max);

f_vector = f_vector(f_index); %frequnecy index 
TW_fft = TW_fft(f_index,:); %

L = length(f_index);


%% (Delay model)
I = scatter_region; 

[sx,sy] = find(I == 1);
% s_x = 105:25:4000;
s_x = sx(1:17)';
% s_x = [102,103,104,105,106,107,108,127,128,129,130,131,132,133,152,153,154,155,156,157,158,177,178,179,180,181,182,183,202,203,204,205,206,207,208,227,228,229,230,231,232,233,252,253,254,255,256,257,258,277,278,279,280,281,282,283];
% s_y = [100, 101,102,103,104,105,106,118,119,120,121,122,123,124,136,137,138,139,140,141,142,154,155,156,157,158,159,160];
s_y = [103];%103,121,139
nx_origin = 95;%s_x(1);
ny_origin =94;
%% Reduce image region
Nx_image=length(s_x);%192;%190;
Ny_image=length(s_y);

M = Nx_image*Ny_image;
DELAY = zeros(L,numelements,M,'single'); 
R = zeros(numelements,M);

x_set = s_x -nx_origin;
y_set = s_y - ny_origin;
[tx,ty] = find(source.p_mask ==1);
%coordinates of tx element
nx_tx = tx - nx_origin;
ny_tx = ty - ny_origin;
%%
m=1;
for nx = x_set %1:Nx_image
    for ny =  y_set %1:Ny_image 
        for j =1:numelements
            Dt = sqrt((nx - nx_tx(j))^2+(ny-ny_tx(j))^2)*dx;
            tof = Dt /c0;
            R(j,m) = Dt;
            DELAY(:,j,m) = exp(-1i*2*pi*f_vector'*tof);  
        end
        m=m+1;
    end
end

 %% generate model matrix for all sensors
A_full = zeros(L*numelements,M,'single');
DELAY_sum = sum(DELAY,2);
for r = 1:numelements
    rec_single = zeros(L,M,'single');
    for t = 1:numelements
        SPR = 1./(4*pi*(R(r,:)));
        new_rec = SPR.*squeeze(DELAY(:,t,:)).* TW_fft(:,t).*squeeze(DELAY(:,r,:));
        rec_single = rec_single + new_rec;
    end
    A_full((r-1)*L+1:L*r,:) = rec_single;
end
