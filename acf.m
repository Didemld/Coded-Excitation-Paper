clc
clear all
close all
%%
load code_sub.mat
aa = Code(80,:);
 aa = [aa, 0,0,0,0];
matched_filter_3= conv(aa,flip(aa));
% plot(matched_filter_3)

%%
load code_sub11.mat
aa = Code(80,:);
 aa = [aa, 0,0,0,0];
matched_filter_4= conv(aa,flip(aa));
% plot(matched_filter_3)



%%
load code_eig.mat
bb = Code(40,:);
bb = [bb, 0,0,0,0];
matched_filter_5 = conv(bb,flip(bb));
% hold on;
% plot(matched_filter_5)

%%
load code_tr.mat
cc = Code(40,:);
cc = [cc, 0,0,0,0];
matched_filter_10 = conv(cc,flip(cc));
% hold on;
% plot(matched_filter_10)
%%
all = [matched_filter_5;matched_filter_10;matched_filter_3;matched_filter_4];
all_norm = all / norm(all);%10*log10(abs(all)./max(all(:)));
xaxis = -5:1:5;
figure;
plot(xaxis,real(all_norm(1,:)),'Linewidth',3);
hold on;
plot(xaxis,real(all_norm(2,:)),'Linewidth',3);
hold on;
plot(xaxis,real(all_norm(3,:)),'Linewidth',3);
hold on;
plot(xaxis,real(all_norm(4,:)),'Linewidth',3);
% axis tight
legend('EIG','TR','SUB','SUB11')
%%
xlabel('Code shift','FontSize', 14)
ylabel('Normalized Amplitude','FontSize', 14)
title('ACF for different bit-length optimized code')