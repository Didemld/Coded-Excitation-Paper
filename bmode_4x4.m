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
close all
%%
tran=5
for i = 1:tran
    load(['I_eig_4x4', num2str(i), '.mat']);
    I_rr4(:,i) = I_matchedfilter;
end 


%%

for i = 1:tran
    load(['I_tr_4x4', num2str(i), '.mat']);
    I_tr4(:,i) = I_matchedfilter;
end 



%%
for i = 1:tran
    load(['I_sa_4x4', num2str(i), '.mat']);
    I_e4(:,i) = I_matchedfilter;
end 

%%
for i = 1:tran
    load(['I_2e11', num2str(i), '.mat']);
    I_be4(:,i) = I_matchedfilter;
end 


%%
%%
Nx_image=295;%190;
Ny_image=120;%72;
nx_origin = 170;
ny_origin =71;%94;


IMAGE_xaxis = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;


rt =20;
t=5;
figure;
subplot(1,t,1)
I_r = mean(I_rr4,2);
Bmode_c= 10*log10(I_r./max(I_r(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_r = im_r(:,11:110);
IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
colorbar
% title('compounded 10 eigenvectors');
title('cvx (Eig-Opt,16 pixels)')
impixelinfo
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
subplot(1,t,2)
I_t = mean(I_tr4,2);
Bmode_c= 10*log10(I_t./max(I_t(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_r = im_r(:,11:110);
IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
colorbar
% title('compounded 10 eigenvectors');
title('cvx (Trace-Opt, 16 pixels)')
impixelinfo
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
subplot(1,t,3);
I_c = mean(I_e4,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_r = im_r(:,11:110);
IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
title('suboptimal (16 pixels)');
colorbar
% title('optimized compounding')
impixelinfo
axis equal; axis tight;

%%
subplot(1,t,4);
I_c = mean(I_be4,2);
Bmode_c= 10*log10(I_c./max(I_c(:)));
im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
im_r = im_r(:,11:110);
IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
caxis([-rt 0])
xlabel('(mm)')
ylabel('(mm)')
title('suboptimal (121 pixels)');
colorbar
% title('optimized compounding')
impixelinfo
axis equal; axis tight;
%%
subplot(1,t,5);
imagesc(IMAGE_yaxis,IMAGE_xaxis,scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin+10:ny_origin+Ny_image-20+10-1));
xlabel('(mm)')
ylabel('(mm)')
colorbar
title('original');
axis equal; axis tight;


%%
figure;
for e = 1:tran
    subplot(1,tran,e)
    Bmode_c= 10*log10(I_rr4(:,e)./max(I_rr4(:,e)));
    im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_r = im_r(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(6-e),'th',' random selection'])
    impixelinfo
    axis equal; axis tight;
end
% 
% %%
% figure;
% for e = 1:tran
%     subplot(1,tran,e)
%     Bmode_e= 10*log10(I_dd(:,e)./max(I_dd(:,e)));
%     imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_e),Ny_image,Nx_image)');
%     caxis([-rt 0])
%     xlabel('(mm)')
%     ylabel('(mm)')
%     title([num2str(e),'th',' eigenvector'])
%     impixelinfo
%     axis equal; axis tight;
% end



%%
figure;
for e = 1:tran
    subplot(1,tran,e)
     Bmode_c= 10*log10(I_tr4(:,e)./max(I_tr4(:,e)));
    im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_r = im_r(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th',' random selection'])
    impixelinfo
    axis equal; axis tight;
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
    Bmode_c= 10*log10(I_e4(:,e)./max(I_e4(:,e)));
    im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';
    im_r = im_r(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);  
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
     Bmode_c= 10*log10(I_be4(:,e)./max(I_be4(:,e)));
    subplot(1,tran,e)
    im_r = reshape(real(Bmode_c),Ny_image,Nx_image)';   
    im_r = im_r(:,11:110);
    IMAGE_yaxis = [(0-(Ny_image-20)/2) ((Ny_image-20)/2-1)]*dx*1e3;
    imagesc(IMAGE_yaxis,IMAGE_xaxis,im_r);
    caxis([-rt 0])
    xlabel('(mm)')
    ylabel('(mm)')
    title([num2str(e),'th',' eigenvector '])
    impixelinfo
    axis equal; axis tight;
end
