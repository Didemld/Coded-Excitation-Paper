%%
%%% This code is to calculate the A_tr matrix based on the selected imaging
%%% The area can be selected at Reduce image region
%%% The code length can be changed by changing code_length
%%% require mat file (initialize_GE, makeG, externel matlab function mtimesx
%%% save A_tr matrix and use this to calcualte the encoding matrix C (theleading eigenvector)

clc
clear all
close all
%%
addpath('/mnt/Data1/didem/didem_tu_delft/k-wave-toolbox-version-1.3/k-Wave')
didem_initialize_k
%% Source = burst
tone_burst_freq = 15*1e6;          % [Hz]
tone_burst_cycles = 3;
sampling_freq = 1/kgrid.dt;     % [Hz]
element_spacing = dx;           % [m]
base_pulse = toneBurst(sampling_freq, tone_burst_freq, tone_burst_cycles);
numelements = num_elements;
%% encoding matrix C
pulse_length = length(base_pulse);
code_length = 2;  %% bit length. could be 3, 5, 7 or even more, in the thesis I chose 5.

%% transform TX to fourier domain

Nt=kgrid.Nt;
TW_zero_padding = zeros(Nt,1);
TW_zero_padding(1:pulse_length) =base_pulse;
tw_fft =fft(TW_zero_padding);
TW_fft = [tw_fft(Nt/2+2:end) ;tw_fft(1:Nt/2+1)];

fs = 1/kgrid.dt; %sampling frequency
f0 = tone_burst_freq;  %transducer frequency
f_vector = fs*(-(Nt/2-1):(Nt/2))/Nt;

f_min = f0 - f0 *.25;
f_max = f0 + f0 *.25;

f_index = find (f_vector >= f_min & f_vector <= f_max);

f_vector = f_vector(f_index); %frequnecy index
f_vector = single(f_vector);
TW_fft = TW_fft(f_index); %

L = length(f_index);
%%
phase = zeros(L,code_length,1);
for i = 1:code_length
    phase(:,i) = exp(-1i*2*pi*f_vector'*(kgrid.dt*(i-1)*pulse_length));
end
%%
p = TW_fft.*phase;
p = single(p);
%% Reduce image region
Nx_image= 102;%389;%260;%204;
Ny_image= 102;%200;%315;%630;%366;%240;%238;



%% (Delay model)
downsample_rate =10; %% further reduce image matrix size
down_x = downsample((1:Nx_image),downsample_rate);
down_y = downsample((1:Ny_image),downsample_rate);
% M = Nx_image*Ny_image;
M = length(down_x)*length(down_y);

DELAY = zeros(L,numelements,M,'single');
R = zeros(numelements,M);
angle = zeros(numelements,M);
nx_origin = 180+(Nx-180)/2-Nx_image/2;%587;%442;%200;%252;%169;200-
ny_origin =Ny/2-Ny_image/2 ;%137;%139;

[tx,ty] = find(source.p_mask ==1);
%coordinates of tx element
nx_tx = tx - nx_origin;
ny_tx = ty - ny_origin;

m=1;
for nx= down_x
    for ny= down_y
        for j =1:numelements
            Dt = sqrt((nx - nx_tx(j))^2+(ny-ny_tx(1 +(j-1)*element_width))^2)*dx;
            tof = Dt /c0;
            angle(j,m) = (abs(nx - nx_tx(j))*dx)/Dt;
            R(j,m) = Dt;
            DELAY(:,j,m) = exp(-1i*2*pi*f_vector'*tof);
        end
        m=m+1;
    end
end

%% bulding transmit green function
g_transmit = permute(DELAY,[3,2,1]);
% g_transmit = reshape(g_transmit,size(g_transmit,1),[],1);
%% buiding receive model
g_receive = permute(DELAY,[2,3,1]);
for j = 1:numelements
    SPR = (angle(j,:))./(4*pi*(R(j,:)));
    g_receive(j,:,:) = SPR'.*squeeze(g_receive(j,:,:));
end
% g_receive = permute(g_receive,[1,3,2]);
% g_recevive = reshape(g_receive,[],size(g_receive,3),1);
c_num = numelements*code_length;

%
addpath('/mnt/Data1/didem/encoding_matrix_code_MacOS/cvx_linux')
cvx_setup;

cvx_solver
cvx_begin
disp('Initializing cvx cost function');


variable lambda;
variable c(c_num,1);
variable C(c_num,c_num);


sumcA=zeros(M);

%% generate model matrix for all sensors
g_transmit = (g_transmit);
g_receive = (g_receive);
p = (p);
A_tr = (zeros(c_num,c_num,'single')); 
for i = 1:c_num
    O_c = (zeros(numelements,code_length,'single'));
    O_c(i) = 1;
    G_1 = makeG(g_receive,g_transmit,O_c,p,L,numelements,M);
    G_1 = permute(G_1,[1,3,2]);
    G_1 = reshape(G_1,[],size(G_1,3),1);
    G_1 = G_1';
    tic
    for j = i:c_num
        O_c = (zeros(numelements,code_length,'single'));
        O_c(j) = 1;
        G_2 = makeG(g_receive,g_transmit,O_c,p,L,numelements,M);

        G_2 = permute(G_2,[1,3,2]);
        G_2 = reshape(G_2,[],size(G_2,3),1);
        A = G_1*G_2; %externel function mtimesx, faster matrix calculation
        A_tr(i,j) = trace(A);
        A_tr(j,i) = A_tr(i,j);
       
        sumcA =sumcA+(C(i,j)+C(j,i))*double(A);
%         A_var(:,:,i,j) = A;
%         A_var(:,:,j,i) = A_var(:,:,i,j);
    end
    toc
    i
end
%% finding the eigenvalues
A_tr = gather(A_tr); %save A_tr to calculate the encoding matrix C.




% Constraints
disp('Initializing cvx constraints');
subject to

% [ C c ; c' 1 ] == semidefinite( size(C,1)+1 );

% c
C == hermitian_semidefinite(length(c));
% Opt constraints

sum(diag(C)) == 1;

% Typical relaxation for discrete problems
sumcA-lambda *eye(M) == hermitian_semidefinite(M);

disp('Finished cvx init');
cvx_end

