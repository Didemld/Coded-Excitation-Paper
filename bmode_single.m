clc
clear all
close all
% %%

addpath('codes_128_3x3_code1')
addpath('/Users/didemdogan/surfdrive/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')

%% create the computational grid
Nx = 512;           % number of grid points in the x (row) direction
Ny = 256;          % number of grid points in the y (column) direction

dx = 0.025e-3;        % grid point spacing in the x direction [m]
dy = 0.025e-3;        % grid point spacing in the y direction [m]


kgrid = kWaveGrid(Nx, dx, Ny, dy);

%% define the properties of the propagation medium

c0 = 1540;                  % [m/s]
rho0 = 1000;                % [kg/m^3]
medium.alpha_coeff = 0.75;  % [dB/(MHz^y cm)]
medium.alpha_power = 1.5;

% medium.alpha_coeff = 0;  % [dB/(MHz^y cm)]
% medium.alpha_power = 0;

scattering_c = 2050;
scattering_rho = 500; 

%% create time array

t_end = (Nx*dx)*2.2/c0;
kgrid.t_array = makeTime(kgrid, c0, [], t_end);

%% transducer array position
num_elements = 128;      % [grid points]
element_width = 1 ; % in grid points
transducer_width = num_elements*element_width; % in grid points  
start_index = Ny/2 - round(transducer_width/2) + 1; 

%% tx source

nx_offset = 25;       
source.p_mask=zeros(Nx,Ny); %same size as the computational grid
%array
source.p_mask(nx_offset,start_index:start_index + transducer_width - 1 )=1;
%% define three scatters
scatter_region = zeros(Nx,Ny);
%%
nx_step = 25;
ny_step=25;
% scatter_region(nx_offset+40:nx_step:end,start_index+5:ny_step:end)=1;

scatter_region(346,128)=1;



% 
sound_speed_map = ones(Nx,Ny).*c0;  % [m/s]
sound_speed_map(scatter_region==1) = scattering_c;
medium.sound_speed = sound_speed_map;

density_map = ones(Nx,Ny).*rho0;
density_map(scatter_region==1) = scattering_rho;
medium.density = density_map;
%% define sensor mask

sensor.mask = zeros(Nx,Ny);
sensor.mask(nx_offset, start_index:start_index + transducer_width - 1) = 1;
% %%
% figure
% imagesc(scatter_region+source.p_mask)

close all

%%
tran=1
for i = 1:tran
    load(['I_eig_s.mat']);
    I_eig(:,i) = I_matchedfilter;

end

for i = 1:tran 
    load(['I_sa_s.mat']);
    I_sa(:,i) = I_matchedfilter;
end 
%%
for i = 1:tran
    load(['I_tra_s1.mat']);
    I_trai(:,i) = I_matchedfilter;
end 

%%
for i = 1:tran
    load(['I_pw_s.mat']);
    I_be(:,i) = I_matchedfilter;
end 


%%
for i = 1:tran
    load(['I_bcs.mat']);
    I_bc(:,i) = I_matchedfilter;
end 

%%
for i = 1:tran
    load(['I_rands.mat']);
    I_rd(:,i) = I_matchedfilter;
end 


%%
Nx_image=295;%190;
Ny_image=120;%72;
nx_origin = 170;
ny_origin =71;%94;

Nx_image=100;%190;
Ny_image=100;%72;
% (Delay model)
nx_origin = 300;
ny_origin =78;%94;


IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;

small_scatter_region = scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin:ny_origin+Ny_image-1);

[ind_row,ind_col] =find(small_scatter_region ==1);
ind_row = ind_row-2:ind_row+2;
ind_col = ind_col-2:ind_col+2;

noise_ind1 = [1:5];
noise_ind2 = [1:5];
%


rt =20;
t=6;
figure;
subplot(2,t,1)
I_r = mean(I_eig,2);
Bmode_c= 10*log10(I_r./max(I_r(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;

imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')


I_mid = reshape(real(I_r./max(I_r(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

% title('compounded 10 eigenvectors');
title('Eig-Opt', ['SNR = ', num2str(SNR)])


% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_r ./ max(I_r(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Eig-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));


axis equal; axis tight;



%%
subplot(2,t,2)
I_t = mean(I_trai,2);
Bmode_c= 10*log10(I_t./max(I_t(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_t./max(I_t(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
% title('compounded 10 eigenvectors');
% title('Trace-Opt', ['SNR = ', num2str(SNR)])
% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_t ./ max(I_t(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

title('Trace-Opt', ['SNR = ', num2str(SNR)])


%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Trace-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));


axis equal; axis tight;


%%
subplot(2,t,3);
I_c = mean(I_sa,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';
I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

title('FIM-Opt', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('FIM-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;


%%
subplot(2,t,4);
I_c = mean(I_bc,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
title('Barker codes', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Barker Codes | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;


%%
subplot(2,t,5);
I_c = mean(I_rd,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
title('Random Codes', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Random Codes | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;



%%
subplot(2,t,6);
I_c = mean(I_be,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(3,3);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
title('Plane-wave', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Plane-wave | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));
axis equal; axis tight;


%%
tran=1
for i = 1:tran
    load(['I_eig_sn.mat']);
    I_eig(:,i) = I_matchedfilter_noise;
end 
% 
% %%
% for i = 1:tran
%     load(['I_2dd', num2str(i), '.mat']);
%     I_dd(:,i) = I_matchedfilter;
% end 


%%

for i = 1:tran 
    load(['I_sa_sn.mat']);
    I_sa(:,i) = I_matchedfilter_noise;
end 

%%

% for i = 1:tran
%     load(['I_trdd', num2str(i), '.mat']);
%     I_trd(:,i) = I_matchedfilter;
% end 

%%
for i = 1:tran
    load(['I_tra_sn1.mat']);
    I_trai(:,i) = I_matchedfilter_noise;
end 

%%
for i = 1:tran
    load(['I_pw_sn.mat']);
    I_be(:,i) = I_matchedfilter_noise;
end 

%%
for i = 1:tran
    load(['I_bcsn.mat']);
    I_bc(:,i) = I_matchedfilter_noise;
end 


%%
for i = 1:tran
    load(['I_randsn.mat']);
    I_rd(:,i) = I_matchedfilter_noise;
end 


%%
Nx_image=295;%190;
Ny_image=120;%72;
nx_origin = 170;
ny_origin =71;%94;

Nx_image=100;%190;
Ny_image=100;%72;
% (Delay model)
nx_origin = 300;
ny_origin =78;%94;


IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;

%%


rt =20;
t=6;
subplot(2,t,7)
I_r = mean(I_eig,2);
Bmode_c= 10*log10(I_r./max(I_r(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_r./max(I_r(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
% title('compounded 10 eigenvectors');
title('Eig-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])
% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_r ./ max(I_r(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Eig-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));
axis equal; axis tight;


%%
subplot(2,t,8)
I_t = mean(I_trai,2);
Bmode_c= 10*log10(I_t./max(I_t(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_t./max(I_t(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
% title('compounded 10 eigenvectors');
title('Trace-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_t ./ max(I_t(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Trace-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;


%%
subplot(2,t,9);
I_c = mean(I_sa,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
title('FIM-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])
% title('optimized compounding')
% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('FIM-Opt | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;



%%
subplot(2,t,10);
I_c = mean(I_bc,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
title('Barker codes/10 dB Noise', ['SNR = ', num2str(SNR)])

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Barker Codes| Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

% title('optimized compounding')

axis equal; axis tight;



%%
subplot(2,t,11);
I_c = mean(I_rd,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

title('Random codes/10 dB Noise', ['SNR = ', num2str(SNR)])
% title('optimized compounding')
% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Random codes | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;


%%
subplot(2,t,12);
I_c = mean(I_be,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
IMAGE_yaxis = [(0-(Ny_image)/2) ((Ny_image)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_mid = reshape(real(I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_mid(ind_row,ind_col);
I_noise = I_mid(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

title('Plane-wave/10 dB Noise', ['SNR = ', num2str(SNR)])
colorbar
% title('optimized compounding')

% ... (your existing plotting code up to SNR calculation)

%% ========================================================================
% Add Resolution Calculation Here
% ========================================================================

% Define the expected scatterer position in the small_scatter_region
% (adjust if your phantom has multiple/moved scatterers)
[target_row, target_col] = find(small_scatter_region == 1);

% Convert to coordinates in the reconstructed image matrix
peak_row = target_row;    % Axial (depth) index
peak_col = target_col;    % Lateral (width) index

% Extract linear-scale image data (normalized)
I_linear = reshape(real(I_c ./ max(I_c(:))), Ny_image, Nx_image)';

% -------------------------------------------------------------------------
% Axial Resolution (FWHM along depth/z-axis)
% -------------------------------------------------------------------------
axial_profile = I_linear(:, peak_col);      % Depth profile through scatterer
[max_val, ~] = max(axial_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = axial_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dx = 0.025 mm/pixel)
axial_res_mm = (last_idx - first_idx) * dx * 1e3; 

% -------------------------------------------------------------------------
% Lateral Resolution (FWHM along width/x-axis)
% -------------------------------------------------------------------------
lateral_profile = I_linear(peak_row, :);    % Width profile through scatterer
[max_val, ~] = max(lateral_profile);
half_max = max_val * 0.5;

% Find FWHM indices
above_threshold = lateral_profile >= half_max;
first_idx = find(above_threshold, 1, 'first');
last_idx = find(above_threshold, 1, 'last');

% Calculate in mm (dy = 0.025 mm/pixel)
lateral_res_mm = (last_idx - first_idx) * dy * 1e3; 

%% ========================================================================
% Update Title with SNR & Resolution
% ========================================================================
title(sprintf('Plane wave | Axial=%.2f mm | Lat=%.2f mm',...
              axial_res_mm, lateral_res_mm));

axis equal; axis tight;




