clc
clear all
close all
% %%

addpath('/mnt/Data1/didem/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/small_array/untitled folder')
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/small_array/buyuk')
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/small_array/untitled 1x1')
addpath('codes_128_3x3_code1')
didem_initialize_k

tran=3;
% ind_noise1 = [1:10];
% ind_noise2= [1:10];

%%
Nx_image=295;%190;
Ny_image=120;%72;
nx_origin = 170;
ny_origin =71;%94;

small_scatter_region = scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin+10:ny_origin+Ny_image-20+10-1);

[ind_row,ind_col] =find(small_scatter_region(:,:) ==1);
% ind_row = [ind_row-3;ind_row-2;ind_row-1;ind_row;ind_row+1;ind_row+2;ind_row+3];
% ind_col = [ind_col-3;ind_col-2;ind_col-1;ind_col;ind_col+1;ind_col+2;ind_col+3];
% close all
[ind_row2,ind_col2] =find(small_scatter_region(1:100,:) ==1);
ind_noise1= ind_row2-13;
ind_noise2=ind_col2-13;
ind_noise1 = [1:5];
ind_noise2 = [1:5];

% ind_noise1 = [ind_noise1-2;ind_noise1-1;ind_noise1;ind_noise1+1;ind_noise1+2];
% ind_noise2 = [ind_noise2-2;ind_noise2-1;ind_noise2;ind_noise2+1;ind_noise2+2];
ind_noise2(ind_noise2<=0) =1;
ind_noise1(ind_noise1<=0) =1;




% ind_noise1(ind_noise1>295) =1;
% ind_noise2(ind_noise2>100) =1;
%%
for i = 1:tran
    load(['I_eig_3x3', num2str(i), '.mat']);
    I_eig_opt(:,i) = I_matchedfilter;
end 

if tran > 5

for i = 6:tran 
    load(['I_eig_3x3a', num2str(i-5), '.mat']);
    I_eig_opt(:,i) = I_matchedfilter;
end 
end
%  I_rr(:,6:10) = I_rr(:,1:5);
% % I_rr(:,7) = I_rr(:,6);
% % I_rr(:,8) = I_rr(:,10);
%%

for i = 1:tran 
    load(['I_tri_3x3', num2str(i), '.mat']);
    I_trace_opt(:,i) = I_matchedfilter;
end 

%%

% for i = 1:tran
%     load(['I_trdd', num2str(i), '.mat']);
%     I_trd(:,i) = I_matchedfilter;
% end 

%%
for i = 1:tran
    load(['I_sa_3x3', num2str(i), '.mat']);
    I_fim_9(:,i) = I_matchedfilter;
end 

%%
for i = 1:tran
    load(['I_sa_11x11', num2str(i), '.mat']);
    I_fim_121(:,i) = I_matchedfilter;
end 





IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;



%%

if tran == 1
    t = 5;
else 
    t = 4;
end


rt =20;
figure;

%%


% ind_noise = [ind_row_noise,ind_col_noise];
IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;


if tran == 1
%    small_scatter_region(ind_noise1,ind_noise2) = 0.5; 
subplot(1,t,5);
imagesc(IMAGE_yaxis,IMAGE_xaxis,small_scatter_region);
xlabel('x (mm)')
ylabel('z (mm)')
colorbar
title('original');
axis equal; axis tight;

end
% sgtitle([num2str(tran), ' transmissions' ])



subplot(1,t,1)
I_r = mean(I_eig_opt,2);


Bmode_c= 10*log10(I_r./max(I_r(:)));
im_eig = reshape(real(Bmode_c),Ny_image,Nx_image)';

im_eig = im_eig(:,11:110);

I_r=I_r./max(I_r(:));

I_mid = reshape((I_r),Ny_image,Nx_image)';

I_mid = I_mid(:,11:110);

% psnr(double(I_mid),small_scatter_region)


I_sig_eig = I_mid(ind_row,ind_col);

I_noise_eig=I_mid(ind_noise1,ind_noise2);

ind_noise_list = [35 9; 9 11;85 17;105 27;129 26; 154 26; 182 27;209 30; 236 30; 263 30;283 31; 37 55; ...
    55 49; 85 53;110 52;136 53; 159 51; 186 52; 209 52; 230 51; 254 52; 276 52;];

I_noise_eig=I_mid(ind_noise_list);

SNR_Eig = 10*log10(mean(mean(abs(I_sig_eig).^2))/mean(mean(abs(I_noise_eig).^2)))


% Calculate the mean intensity of signal and noise
mu_signal = mean(mean(abs(I_sig_eig)));    % Mean signal intensity
mu_noise = mean(mean(abs(I_noise_eig)));  % Mean noise intensity

% Calculate the standard deviation of noise
sigma_noise = std(abs(I_noise_eig(:)));   % Standard deviation of noise intensities

% Compute CNR
CNR_Eig = (mu_signal - mu_noise) / sigma_noise;

% im_eig(ind_noise1, ind_noise2) =0;

% imagesc(IMAGE_yaxis,IMAGE_xaxis,im_eig);
imagesc(im_eig);


caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
%colorbar
% title('compounded 10 eigenvectors');
title('Eig-Opt, 9 pixels', ['SNR = ', num2str(SNR_Eig)])
axis equal; axis tight;


%%
subplot(1,t,2)
I_t = mean(I_trace_opt,2);


Bmode_c= 10*log10(I_t./max(I_t(:)));
im_tra = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_tra = im_tra(:,11:110);

I_t=I_t./max(I_t(:));
I_mid = reshape((I_t),Ny_image,Nx_image)';

I_mid = I_mid(:,11:110);
% psnr(double(I_mid),small_scatter_region)


I_sig_tra = I_mid(ind_row,ind_col);
I_noise_tra=I_mid(ind_noise1,ind_noise2);

ind_noise_list = [8 8; 34 8; 60 8; 106 8; 130 3; 154 4; 32 23; 54 22; 104 28; 135 27; 156 25; 175 26; ...
    208 26;238 27; 255 29;254 29; 284 29; 33 46;58 47;79 49; 108 49; 138 52; 156 52; 180 47; 195 50;];

I_noise_tra=I_mid(ind_noise_list);


SNR_Tra = 10*log10(mean(mean(abs(I_sig_tra).^2))/mean(mean(abs(I_noise_tra).^2)))

% Calculate the mean intensity of signal and noise
mu_signal = mean(mean(abs(I_sig_tra)));    % Mean signal intensity
mu_noise = mean(mean(abs(I_noise_tra)));  % Mean noise intensity

% Calculate the standard deviation of noise
sigma_noise = std(abs(I_noise_tra(:)));   % Standard deviation of noise intensities

% Compute CNR
CNR_Tra = (mu_signal - mu_noise) / sigma_noise;


IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
% imagesc(IMAGE_yaxis,IMAGE_xaxis,im_tra);
imagesc(im_tra);

caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
%colorbar
% title('compounded 10 eigenvectors');
title('Trace-Opt, 9 pixels',['SNR = ', num2str(SNR_Tra)])
axis equal; axis tight;


%%
subplot(1,t,3);
I_c = mean(I_fim_9,2);


Bmode_c= 10*log10(I_c./max(I_c(:)));
im_fim9 = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_fim9= im_fim9(:,11:110);

I_c=I_c./max(I_c(:));

I_mid = reshape((I_c),Ny_image,Nx_image)';

I_mid = I_mid(:,11:110);

% psnr(double(I_mid),small_scatter_region)

I_sig_fim9 = I_mid(ind_row,ind_col);

I_noise_fim9=I_mid(ind_noise1,ind_noise2);

ind_noise_list = [3 8; 31 6; 56 9; 76 7; 100 10; 132 6; 158 9; 175 9; 204 11;238 2; 7 32; 33 32; ...
    55 30; 56 31;79 32; 102 29; 130 25; 159 25; 185 24;212 25;244 25;8 46;32 47;55 45;78 49; 106 49; 129 49;...
    158 51;209 44; 232 45];
I_noise_fim9=I_mid(ind_noise_list);


SNR_Fim9 = 10*log10(mean(mean(abs(I_sig_fim9).^2))/mean(mean(abs(I_noise_fim9).^2)))

% Calculate the mean intensity of signal and noise
mu_signal = mean(mean(abs(I_sig_fim9)));    % Mean signal intensity
mu_noise = mean(mean(abs(I_noise_fim9)));  % Mean noise intensity

% Calculate the standard deviation of noise
sigma_noise = std(abs(I_noise_fim9(:)));   % Standard deviation of noise intensities

% Compute CNR
CNR_Fim9 = (mu_signal - mu_noise) / sigma_noise;



IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
% imagesc(IMAGE_yaxis,IMAGE_xaxis,im_fim9);
imagesc(im_fim9);
caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
title('FIM-Opt, 9 pixels',['SNR = ', num2str(SNR_Fim9)]);
%colorbar
% title('optimized compounding')
axis equal; axis tight;

%%
subplot(1,t,4);
I_c = mean(I_fim_121,2);



Bmode_c= 10*log10(I_c./max(I_c(:)));
im_fim121 = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_fim121 = im_fim121(:,11:110);

I_c=I_c./max(I_c(:));

I_mid = reshape((I_c),Ny_image,Nx_image)';

I_mid = I_mid(:,11:110);

% psnr(double(I_mid),small_scatter_region)

I_sig_fim121 = I_mid(ind_row,ind_col);
I_noise_fim121=I_mid(ind_noise1,ind_noise2);


ind_noise_list= [7 7 ;29 7;54 6; 78 6; 104 7;129 8; 154 7; 175 6; 199 5; 202 4; 226 3; 9 23; 34 22;57 22;...
    58 47; 82 49; 108 48; 132 48; 156 48;];

I_noise_fim121=I_mid(ind_noise_list);



SNR_Fim121 = 10*log10(mean(mean(abs(I_sig_fim121).^2))/mean(mean(abs(I_noise_fim121).^2)))

% Calculate the mean intensity of signal and noise
mu_signal = mean(mean(abs(I_sig_fim121)));    % Mean signal intensity
mu_noise = mean(mean(abs(I_noise_fim121)));  % Mean noise intensity

% Calculate the standard deviation of noise
sigma_noise = std(abs(I_noise_fim121(:)));   % Standard deviation of noise intensities

% Compute CNR
CNR_Fim121 = (mu_signal - mu_noise) / sigma_noise;


IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
% imagesc(IMAGE_yaxis,IMAGE_xaxis,im_fim121);
imagesc(im_fim121);

caxis([-rt 0])
xlabel('x (mm)')
ylabel('z (mm)')
title('FIM-Opt, 121 pixels',['SNR = ', num2str(SNR_Fim121)]);
colorbar
% title('optimized compounding')
axis equal; axis tight;


%%
  
% % 
% % %%
% % figure;
% % for e = 1:tran
% %     subplot(1,tran,e)
% %     Bmode_e= 10*log10(I_dd(:,e)./max(I_dd(:,e)));
% %     imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
% %     caxis([-rt 0])
% %     xlabel('(mm)')
% %     ylabel('(mm)')
% %     title([num2str(e),'th',' eigenvector'])
% %     impixelinfo
% %     axis equal; axis tight;
% % end
% 
% 
% 
%%
if tran == 5
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_c= 10*log10(I_eig_opt(:,e)./max(I_eig_opt(:,e)));
    im_eig = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_eig = im_eig(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_eig);
    caxis([-rt 0])
    xlabel('x (mm)')
    ylabel('z (mm)')
    title([num2str(e),'th',' random selection'])
    axis equal; axis tight;
    if e == tran
        colorbar
    end
    
end  

figure;
for e = 1:tran
    subplot(1,tran,e)
     Bmode_c= 10*log10(I_trace_opt(:,e)./max(I_trace_opt(:,e)));
    im_eig = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_eig = im_eig(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_eig);
    caxis([-rt 0])
    xlabel('x (mm)')
    ylabel('z (mm)')
    title([num2str(e),'th',' random selection'])
    axis equal; axis tight;
    if e == tran
        colorbar
    end
end

% %%
% figure;
% for e = 1:tran
%     subplot(1,tran,e)
%     Bmode_e= 10*log10(I_trd(:,e)./max(I_trd(:,e)));
%     imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
%     caxis([-rt 0])
%     xlabel('(mm)')
%     ylabel('(mm)')
%     title([num2str(e),'th',' eigenvector '])
%     impixelinfo
%     axis equal; axis tight;
% end


%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_c= 10*log10(I_fim_9(:,e)./max(I_fim_9(:,e)));
    im_eig = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_eig = im_eig(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_eig);  
    caxis([-rt 0])
    xlabel('x (mm)')
    ylabel('z (mm)')
    title([num2str(e),'th',' eigenvector '])
    axis equal; axis tight;
    
        if e == tran
        colorbar
    end
end
%%
figure;
for e = 1:tran
     Bmode_c= 10*log10(I_fim_121(:,e)./max(I_fim_121(:,e)));
    subplot(1,tran,e)
    im_eig = reshape(real(Bmode_c),Ny_image,Nx_image)';   
    im_eig = im_eig(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_eig);
    caxis([-rt 0])
    xlabel('x (mm)')
    ylabel('z (mm)')
    title([num2str(e),'th',' eigenvector '])
    axis equal; axis tight;
        if e == tran
        colorbar
    end
end

end
impixelinfo

%%%%%%%%

% xaxis = linspace(IMAGE_xaxis(1),IMAGE_xaxis(2),295);
% colu=40;
% 
% figure 
% plot(xaxis,im_eig(:,colu),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
% 
%     
% plot(xaxis,im_tra(:,colu),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
%    
% 
% plot(xaxis,im_fim9(:,colu),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
%  
%    
%    plot(xaxis,im_fim121(:,colu),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
% 
%    
%       
%    plot(xaxis,small_scatter_region(:,colu),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
%     legend('Eig-Opt','Tra-Opt', 'FIM-Opt9', 'FIM-Opt121')
%     
%     datacursormode on
   
% %%
% Bmode_rn= 10*log10(I_rn./max(I_rn(:)));
% im_rn = reshape(real(Bmode_rn),Ny_new,Nx_new)';
% imagesc(im_rn);
% figure;
% plot(xaxis,im_rn(:,46),'LineWidth',2);hold on;
%     xlabel('z (mm)','FontSize', 18)
%     ylabel('dB','FontSize', 18)
%    axis equal; axis tight;
% nexttile
% imagesc(IMAGE_yaxis2,IMAGE_xaxis2,scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin:ny_origin+Ny_image-1));
% axis equal; axis tight;

% rt = 40
% figure
% subplot(1,3,1)
% imagesc(real(im_3n))
% caxis([-rt 0])
% 
% subplot(1,3,2)
% imagesc(real(im_bn))
% caxis([-rt 0])
% 
% subplot(1,3,3)
% imagesc(real(im_on))
% caxis([-rt 0])




