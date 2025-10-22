%% create the computational grid
Nx = 400;           % number of grid points in the x (row) direction
Ny = 64;          % number of grid points in the y (column) direction

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
num_elements = 32;      % [grid points]
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
nx_step = 50;
ny_step=50;
scatter_region(nx_offset+40:nx_step:end,num_elements)=1;



%%
% scatter_region = zeros(Nx,Ny);
% scatter_region(Nx/2, Ny/2) = 1;  %(128,128)
% scatter_region(Nx/2+20, Ny/2) = 1; %(128,148)
% scatter_region(Nx/2-20, Ny/2) = 1;
% scatter_region(Nx/2, Ny/2+20) = 1; %(128,148)
% scatter_region(Nx/2, Ny/2-20) = 1; %(128,108)
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
%%
figure
imagesc(scatter_region+source.p_mask)
impixelinfo