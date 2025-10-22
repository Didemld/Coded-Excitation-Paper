clc
clear all
close all

%%
load C_opt.mat
bb = Code(80,:)/sqrt(10);
% bb=bb/norm(bb);
matched_filter_5 = conv(bb,flip(bb));
% hold on;
% plot(matched_filter_5)

%%
cc = [1 1 1 -1 1];
matched_filter_bar = conv(cc,flip(cc));
% hold on;
% plot(matched_filter_10)
%%
all_norm= [matched_filter_5;matched_filter_bar];
xaxis = -4:1:4;
all_norm = all_norm/norm(all_norm);
figure;
plot(xaxis,real(all_norm(1,:)),'Linewidth',3);
hold on;
plot(xaxis,real(all_norm(2,:)),'Linewidth',3);
% axis tight

%%
xlabel('Code shift','FontSize', 14)
ylabel('Normalized Amplitude','FontSize', 14)
title('ACF for optimized code and Barker code','FontSize', 14)
legend('Optimized Code', 'Barker Code','FontSize', 11)