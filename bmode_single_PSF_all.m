clc
clear all
close all
% %%

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


I_mide = reshape((I_r./max(I_r(:))),Ny_image,Nx_image)';

I_sig= I_mide(ind_row,ind_col);
I_noise = I_mide(noise_ind1,noise_ind2);
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





% subplot(1,t,2)
% I_d = mean(I_dd,2);
% Bmode_c= 10*log10(I_d./max(I_d(:)));
% imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
% caxdddis([-rt 0])
% xlabel('(mm)')
% ylabel('(mm)')xx
% colorbar
% % title('compounded 10 eigenvectors');
% title('cvx-E (rank-5, 9 pixels)')
% impixelinfo
% axis equal; axis tight;

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
I_midt = reshape((I_t./max(I_t(:))),Ny_image,Nx_image)';

I_sig= I_midt(ind_row,ind_col);
I_noise = I_midt(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)


% title('compounded 10 eigenvectors');
title('Trace-Opt', ['SNR = ', num2str(SNR)])

axis equal; axis tight;


% %%
% subplot(1,t,4)
% I_t = mean(I_trd,2);
% Bmode_c= 10*log10(I_t./max(I_t(:)));
% imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
% caxis([-rt 0])
% xlabel('(mm)')
% ylabel('(mm)')
% colorbar
% % title('compounded 10 eigenvectors');
% title('cvx-A (rank-5, 9 pixels)')
% impixelinfo
% axis equal; axis tight;
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
I_midf = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';
I_sig= I_midf(ind_row,ind_col);
I_noise = I_midf(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise

CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

title('FIM-Opt', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

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
I_midb = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midb(ind_row,ind_col);
I_noise = I_midb(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)


title('Barker codes', ['SNR = ', num2str(SNR)])
% title('optimized compounding')
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
I_midr = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midr(ind_row,ind_col);
I_noise = I_midr(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

title('Random Codes', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

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
I_midp = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midp(ind_row,ind_col);
I_noise = I_midp(3,3);


SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)


title('Plane-wave', ['SNR = ', num2str(SNR)])
% title('optimized compounding')
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
I_miden = reshape((I_r./max(I_r(:))),Ny_image,Nx_image)';

I_sig= I_miden(ind_row,ind_col);
I_noise = I_miden(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise

CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

% title('compounded 10 eigenvectors');
title('Eig-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])

axis equal; axis tight;



% subplot(1,t,2)
% I_d = mean(I_dd,2);
% Bmode_c= 10*log10(I_d./max(I_d(:)));
% imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
% caxis([-rt 0])
% xlabel('(mm)')
% ylabel('(mm)')
% colorbar
% % title('compounded 10 eigenvectors');
% title('cvx-E (rank-5, 9 pixels)')
% impixelinfo
% axis equal; axis tight;

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
I_midtn = reshape((I_t./max(I_t(:))),Ny_image,Nx_image)';

I_sig= I_midtn(ind_row,ind_col);
I_noise = I_midtn(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)


% title('compounded 10 eigenvectors');
title('Trace-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])

axis equal; axis tight;



% %%
% subplot(1,t,4)
% I_t = mean(I_trd,2);
% Bmode_c= 10*log10(I_t./max(I_t(:)));
% imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
% caxis([-rt 0])
% xlabel('(mm)')
% ylabel('(mm)')
% colorbar
% % title('compounded 10 eigenvectors');
% title('cvx-A (rank-5, 9 pixels)')
% impixelinfo
% axis equal; axis tight;
%%
subplot(2,t,9);
I_c = mean(I_sa,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
I_midfn = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midfn(ind_row,ind_col);
I_noise = I_midfn(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise

CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

title('FIM-Opt/10 dB Noise', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

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
I_midbn = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midbn(ind_row,ind_col);
I_noise = I_midbn(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise

CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

title('Barker codes/10 dB Noise', ['SNR = ', num2str(SNR)])

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
I_midrn = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midrn(ind_row,ind_col);
I_noise = I_midrn(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)


title('Random codes/10 dB Noise', ['SNR = ', num2str(SNR)])
% title('optimized compounding')

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
I_midpn = reshape((I_c./max(I_c(:))),Ny_image,Nx_image)';

I_sig= I_midpn(ind_row,ind_col);
I_noise = I_midpn(noise_ind1,noise_ind2);
SNR = 10*log10(mean(mean(abs(I_sig).^2))/mean(mean(abs(I_noise).^2)))
CNR_other = abs(mu_sig-mu_noise)/ sqrt(sigma_sig^2+sigma_noise^2)

mu_sig = mean(I_sig(:));
mu_noise = mean(I_noise(:));
sigma_sig = std(I_sig(:));
sigma_noise = std(I_noise(:));

CNR = abs(mu_sig - mu_noise) / sigma_noise

title('Plane-wave/10 dB Noise', ['SNR = ', num2str(SNR)])
colorbar
% title('optimized compounding')

axis equal; axis tight;


%%PSFs

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_mide);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Eig'
FWHM_L
FWHM_A


% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midt);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Tra'
FWHM_L
FWHM_A


% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midf);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'FIM'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midb);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Barker'
FWHM_L
FWHM_A


% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midr);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Random'
FWHM_L
FWHM_A


% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midp);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'PW'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_miden);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'NOISY'

'Eig'
FWHM_L
FWHM_A


% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midtn);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');
'Tra'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midfn);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'FIM'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midbn);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Barker'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midrn);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'Random'
FWHM_L
FWHM_A

% Calculate Axial and Lateral PSF using the updated function
[lateralPSF, axialPSF, FWHM_L, FWHM_A] = CalPSF_Auto(I_midpn);

% Plot the results
figure;
subplot(1, 2, 1);
plot(lateralPSF);
title(['Lateral PSF, FWHM = ', num2str(FWHM_L), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');


subplot(1, 2, 2);
plot(axialPSF);
title(['Axial PSF, FWHM = ', num2str(FWHM_A), ' [mm]']);
xlabel('Index');
ylabel('Intensity [dB]');

'PW'
FWHM_L
FWHM_A