clc
clear all
close all
% %%

addpath('/mnt/Data1/didem/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/small_array/untitled folder')
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/small_array/untitled 1x1')
addpath('buyuk')
addpath('codes_128_3x3_code1')
didem_2k
%%
tran=5
for i = 1:tran
    load(['I_3x3_EIG', num2str(i), '.mat']);
    I_rr(:,i) = I_matchedfilter_noise;
end

%%

for i = 1:tran
    load(['I_3x3_TR', num2str(i), '.mat']);
    I_dd(:,i) = I_matchedfilter;
end

%%
for i = 1:tran
    load(['I_3x3e', num2str(i), '.mat']);
    I_e(:,i) = I_matchedfilter;
end

%%
for i = 1:tran
    load(['I_3x3eB', num2str(i), '.mat']);
    I_be(:,i) = I_matchedfilter;
end


%% Reduce image region
Nx_image=100;%190;
Ny_image=70;%72;
%% (Delay model)
nx_origin = 150;
ny_origin =29;%94;

IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;

%%
rt =25;
t=5;
figure;
subplot(1,t,1)
I_r = mean(I_rr,2);
Bmode_c= 10*log10(I_r./max(I_r(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
colorbar
% title('compounded 10 eigenvectors');
title('cvx optimized-random rounding')
impixelinfo
axis equal; axis tight;

subplot(1,t,2)
I_d = mean(I_dd,2);
Bmode_c= 10*log10(I_d./max(I_d(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_c),Ny_image,Nx_image)');
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
colorbar
% title('compounded 10 eigenvectors');
title('cvx optimized-rank 1')
impixelinfo
axis equal; axis tight;
%%
subplot(1,t,3);
I_c = mean(I_e,2);
Bmode_cn= 10*log10(I_c./max(I_c(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_cn),Ny_image,Nx_image)');
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
title('suboptimal-9');
colorbar
% title('optimized compounding')
impixelinfo
axis equal; axis tight;

%%
subplot(1,t,4);
I_c = mean(I_be,2);
Bmode_cn= 10*log10(I_c./max(I_c(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_cn),Ny_image,Nx_image)');
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
title('suboptimal-121');
colorbar
% title('optimized compounding')
impixelinfo
axis equal; axis tight;
%%
subplot(1,t,5);
imagesc(IMAGE_yaxis,IMAGE_xaxis,scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin:ny_origin+Ny_image-1));
xlabel('(mm)')
ylabel('(mm)')
colorbar
title('original');
axis equal; axis tight;


%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_e= 10*log10(I_rr(:,e)./max(I_rr(:,e)));
    imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th','random selection'])
    impixelinfo
    axis equal; axis tight;
end


%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_e= 10*log10(I_dd(:,e)./max(I_dd(:,e)));
    imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th',' eigenvector '])
    impixelinfo
    axis equal; axis tight;
end


%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_e= 10*log10(I_e(:,e)./max(I_e(:,e)));
    imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th',' eigenvector '])
    impixelinfo
    axis equal; axis tight;
end
%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_e= 10*log10(I_be(:,e)./max(I_be(:,e)));
    imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th',' eigenvector '])
    impixelinfo
    axis equal; axis tight;
end
