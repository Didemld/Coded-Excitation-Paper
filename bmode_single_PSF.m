clc
clear all
close all

% Add necessary paths
addpath('codes_128_3x3_code1')
addpath('/Users/didemdogan/surfdrive/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')

didem_initialize_single
close all

%%
tran=1
for i = 1:tran
    load(['I_eig_s.mat']);
    I_eig(:,i) = I_matchedfilter;
end 
% 
% %%
% for i = 1:tran
%     load(['I_2dd', num2str(i), '.mat']);
%     I_dd(:,i) = I_matchedfilter;
% end 


%%

for i = 1:tran 
    load(['I_sa_s.mat']);
    I_sa(:,i) = I_matchedfilter;
end 

%%

% for i = 1:tran
%     load(['I_trdd', num2str(i), '.mat']);
%     I_trd(:,i) = I_matchedfilter;
% end 

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
SNR2 = 10*log10(mean(mean(I_sig))^2/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR2 = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)
% Calculate CNR
% CNR_other = abs(mu_sig - mu_noise) / sqrt(sigma_sig^2 + sigma_noise^2)


% Define the position for the text (slightly offset from bottom right corner)
xPos = Nx_image- 20; % 10 pixels offset from the right
yPos = Ny_image - 20; % 10 pixels offset from the bottom

% Write the value to the image
text(xPos, yPos, ['SNR = ', num2str(SNR)], 'Color', 'white', 'FontSize', 12, 'HorizontalAlignment', 'right', 'VerticalAlignment', 'bottom');


% title('compounded 10 eigenvectors');
title('Eig-Opt', ['SNR = ', num2str(SNR)])

axis equal; axis tight;







% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_mid);

% Plot the results
figure;
subplot(2, 2, 1);
imagesc(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [units]']);
xlabel('Index');
ylabel('Intensity [dB]');
colormap gray;
colorbar;

subplot(2, 2, 2);
imagesc(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [units]']);
xlabel('Index');
ylabel('Intensity [dB]');
colormap gray;
colorbar;