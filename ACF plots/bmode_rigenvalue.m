clc
clear all
close all
%%
addpath('C:\Users\lzhu\Documents\MATLAB\GE probe')
initialize_GE
%%
load('I_3n.mat');
load('I_3.mat');
load('I_5n.mat');
load('I_5.mat');
load('I_10n.mat');
load('I_10.mat');

%%
f0 = 2.8*1e6;  
dx_new = (c0/f0)/2;
dy_new = dx_new;
Nx_image=931;
Ny_image=275;
Nx_new = floor((Nx_image - 1) * (dx / dx_new));
Ny_new = floor((Ny_image - 1) * (dy / dy_new));
nx_origin = 86;
ny_origin =116;

IMAGE_xaxis = [(nx_origin-nx_offset)*dx*1e3 (nx_origin-nx_offset)*dx*1e3+Nx_new*dx_new*1e3];
IMAGE_yaxis = [(0-Ny_new/2) (Ny_new/2-1)]*dx_new*1e3;

IMAGE_xaxis2 = [(nx_origin-nx_offset) (nx_origin+Nx_image-nx_offset)]*dx*1e3;
IMAGE_yaxis2 = [(0-Ny_image/2) (Ny_image/2-1)]*dx*1e3;
%% image with noise
t = tiledlayout(1,3,'TileSpacing','Compact');
t.Padding = 'compact';
nexttile
Bmode_3n= 10*log10(I_3n./max(I_3n(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_3n),Ny_new,Nx_new)');
caxis([-40 0])
xlabel('x (mm)','FontSize', 12)
ylabel('z (mm)','FontSize', 12)
title('3-bit optimized codes')
axis equal; axis tight;

nexttile
Bmode_5n= 10*log10(I_5n./max(I_5n(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_5n),Ny_new,Nx_new)');
caxis([-40 0])
xlabel('x (mm)','FontSize', 12)
ylabel('z (mm)','FontSize', 12)
title('5-bit optimized codes')
axis equal; axis tight;

nexttile
Bmode_10n= 10*log10(I_10n./max(I_10n(:)));
imagesc(IMAGE_yaxis,IMAGE_xaxis,reshape(real(Bmode_10n),Ny_new,Nx_new)');
caxis([-40 0])
xlabel('x (mm)','FontSize', 12)
ylabel('z (mm)','FontSize', 12)
title('10-bit optimized codes')
axis equal; axis tight;

% nexttile
% imagesc(IMAGE_yaxis2,IMAGE_xaxis2,scatter_region(nx_origin:nx_origin+Nx_image-1,ny_origin:ny_origin+Ny_image-1));
% axis equal; axis tight;
