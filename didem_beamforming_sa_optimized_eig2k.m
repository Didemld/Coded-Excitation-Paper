clc
clear all
close all
%%
addpath('/mnt/Data1/didem/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')
didem_2k
%% Source = burst
tone_burst_freq = 15*1e6;          % [Hz]
tone_burst_cycles = 3;
sampling_freq = 1/kgrid.dt;     % [Hz]
element_spacing = dx;           % [m]
base_pulse = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);
numelements = num_elements;
%% encoding matrix C
pulse_length = length(base_pulse);
code_length =3;
load A_tr_192_3x3_code3.mat
load A_var_192_3x3_code3.mat
A_real = real(A_tr);
[V,D] = eig(A_real);  % calculate eigenvectors and eigenvalues of A
[lambda, idx] = max(diag(D));  % find largest eigenvalue and its index
for tran = 1:5
    didem_2k
    load C_192_3x3_code3.mat
    [c_best,C_best,trinv] = C_random_rounding_eig(C,A_var,A_tr,100);
    
    
    c_real = c_best;
    % c_real = V(:, idx-tran+1);  % extract corresponding complex eigenvector
    
    
    c_normalized = c_real / norm(c_real);  % normalize to satisfy constraint
    C = reshape(c_normalized,numelements,code_length);
    Code = C*100;
    %%
    excitation_length = length(base_pulse);
    excitation_signals = zeros(excitation_length*code_length,transducer_width);
    for m = 1:numelements
        for k = 1:code_length
            excitation_signals((k-1)*excitation_length+1:k*excitation_length,(m-1)*element_width+1:m*element_width) = repmat((base_pulse*Code(m,k))',1,element_width);
        end
    end
    excitation_signals = excitation_signals.*10;
    
    %%
    filtered_signals=zeros(length(kgrid.t_array),transducer_width);
    for i=1:transducer_width
        excitation_signal_zeropad = zeros(length(kgrid.t_array),1);
        excitation_signal_zeropad(1:size(excitation_signals,1))=excitation_signals(:,i);
        filtered_signals(:,i) = filterTimeSeries(kgrid, medium, excitation_signal_zeropad);
    end
    %% Source = burst
    source.p = (filtered_signals)';
    %% RF with scatter
    % define the acoustic parameters to record
    sensor.record = {'p'};
    % run the simulation
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
    RF_upsampled_scatter = sensor_data.p;
    
    %% RF without scattering
    scatter_region = zeros(Nx,Ny);
    sound_speed_map = ones(Nx,Ny).*c0;  % [m/s]
    sound_speed_map(scatter_region==1) = scattering_c;
    medium.sound_speed = sound_speed_map;
    
    density_map = ones(Nx,Ny).*rho0;
    density_map(scatter_region==1) = scattering_rho;
    medium.density = density_map;
    
    sensor_data = kspaceFirstOrder2D(kgrid, medium, source, sensor);
    
    RF_upsampled_noScatter = sensor_data.p;
    
    %% downsample RF
    RF_noScatter =zeros(num_elements,kgrid.Nt);
    RF_scatter = zeros(num_elements,kgrid.Nt);
    for i=1:num_elements
        RF_noScatter(i,:) = sum(RF_upsampled_noScatter((i-1)*element_width+1:i*element_width,:),1);
        RF_scatter(i,:) = sum(RF_upsampled_scatter((i-1)*element_width+1:i*element_width,:),1);
    end
    %% transform TX to fourier domain
    
    Nt=kgrid.Nt-2;
    TW_zero_padding = zeros(Nt,numelements);
    tw_fft = zeros(Nt,1);
    TW_fft = zeros(Nt,numelements);
    
    for i = 1:numelements
        TW_zero_padding(:,i) = source.p((i-1)*element_width+1,1:Nt);
        tw_fft =fft(TW_zero_padding(:,i));
        TW_fft(:,i) = [tw_fft(Nt/2+2:end) ;tw_fft(1:Nt/2+1)];
    end
    
    fs = 1/kgrid.dt; %sampling frequency
    f0 = tone_burst_freq;  %transducer frequency
    f_vector = fs*(-(Nt/2-1):(Nt/2))/Nt;
    
    f_min = f0 - f0 *.3;
    f_max = f0 + f0 *.3;
    
    f_index = find (f_vector >= f_min & f_vector <= f_max);
    
    f_vector = f_vector(f_index); %frequnecy index
    TW_fft = TW_fft(f_index,:); %
    
    L = length(f_index);
    
    %% Reduce image region
    Nx_image=100;%190;
    Ny_image=70;%72;
    %% (Delay model)
    nx_origin = 150;
    ny_origin =29;%94;
    M = Nx_image*Ny_image;
    
    DELAY = zeros(L,numelements,M,'single');
    R = zeros(numelements,M);
    angle = zeros(numelements,M);
    
    [tx,ty] = find(source.p_mask ==1);
    %coordinates of tx element
    nx_tx = tx - nx_origin;
    ny_tx = ty - ny_origin;
    m=1;
    for nx=1:Nx_image
        for ny=1:Ny_image
            for j =1:numelements
                Dt = sqrt((nx - nx_tx(j))^2+(ny-ny_tx(j))^2)*dx;
                tof = Dt /c0;
                angle(j,m) = (abs(nx - nx_tx(j))*dx)/Dt;
                R(j,m) = Dt;
                DELAY(:,j,m) = exp(-1i*2*pi*f_vector'*tof);
            end
            m=m+1;
        end
    end
    
    %% generate model matrix for all sensors
    A_full = zeros(L*numelements,M,'single');
    for r = 1:numelements
        SPR = (angle(j,:))./(4*pi*(R(j,:)));
        temp = bsxfun(@times,DELAY, reshape(SPR,1,1,M));
        temp = bsxfun(@times, temp,TW_fft);
        A_full((r-1)*L+1:L*r,:) = squeeze(sum(bsxfun(@times, temp,DELAY(:,r,:)),2));
    end
    
    %% add TGC
    % create radius variable assuming that t0 corresponds to the middle of the
    % input signal
    t0 = pulse_length * kgrid.dt*code_length / 2;
    r = c0 * ( (1:length(kgrid.t_array)) * kgrid.dt/ 2 - t0) ;    % [m]
    
    % define absorption value and convert to correct units
    tgc_alpha_db_cm = medium.alpha_coeff * (tone_burst_freq * 1e-6)^medium.alpha_power;
    tgc_alpha_np_m = tgc_alpha_db_cm / 8.686 * 100;
    
    % create time gain compensation function based on attenuation value and
    % round trip distance
    tgc = exp(tgc_alpha_np_m *2 * r);
    
    %% Receive model R
    
    RF_mute = RF_scatter-RF_noScatter;
    %apply the time gain compensation to each of the scan lines
    RF_tgc = bsxfun(@times, tgc, RF_mute);
    
    %%
    rf = zeros(L,numelements);
    for i = 1:numelements
        rf_fft = transpose(fft(RF_tgc(i,1:Nt)));
        rf_fft =  [rf_fft(Nt/2+2:end) ;rf_fft(1:Nt/2+1)];
        rf(:,i) =  rf_fft(f_index);
    end
    %% add noise
    %     n1 = wgn(length(tgc),1,-4)/400;
    load noise_10db.mat
    
    for v = 1:numelements
        RF_noise(v,:)= RF_mute(v,:); %+n1';
    end
    RF_noisetgc = bsxfun(@times, tgc, RF_noise);
    rf_noise = zeros(L,numelements);
    for i = 1:numelements
        rf_fft_noise = transpose(fft(RF_noisetgc(i,1:Nt)));
        rf_fft_noise =  [rf_fft_noise(Nt/2+2:end) ;rf_fft_noise(1:Nt/2+1)];
        rf_noise(:,i) =  rf_fft_noise(f_index);
    end
    %% matched filter after adding noise
    I_matchedfilter_noise = A_full'*reshape(rf_noise,[],1);
    
    %% Matched filtering
    I_matchedfilter = A_full'*reshape(rf,[],1);
    filename1 =['I_3x3_EIG' num2str(tran)];
    save([filename1 '.mat'],'I_matchedfilter');
    filename2 =['I_3x3_EIGn' num2str(tran)];
    save([filename2 '.mat'],'I_matchedfilter_noise');
end
%%
rt =18;
didem_2k
IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;
figure;
subplot(1,3,1);
Bmode = 10*log10(I_matchedfilter./max(I_matchedfilter(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode),Ny_image,Nx_image)');
caxis([-rt 0])
colorbar
title('optimized')
xlabel('(mm)')
ylabel('(mm)')
impixelinfo
axis equal; axis tight;
subplot(1,3,2);
Bmode_noise = 10*log10(I_matchedfilter_noise./max(I_matchedfilter_noise(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_noise),Ny_image,Nx_image)');
caxis([-rt 0])
title('SNR = 10dB')
impixelinfo
axis equal; axis tight;
colorbar
subplot(1,3,3)
imagesc(IMAGE_yaxis,IMAGE_xaxis,scatter_region(nx_origin:nx_origin+Nx_image,ny_origin:ny_origin+Ny_image));
axis equal; axis tight;
colorbar