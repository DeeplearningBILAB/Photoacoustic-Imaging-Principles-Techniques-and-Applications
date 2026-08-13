% Generating simulated PACT data (Sinogram) using k-wave
% Full clockwise circular scanning with UST
% Point target numerical phantom

% author: Manojit Pramanik
% date: 05 August 2026
% k-Wave Version 1.4.1 and Matlab R2026a Update 4 (26.1.0.3312084)

clear all; close all; clc;
%% parameters
medium.sound_speed = 1500; % Speed of sound in water [m/s]
outputfilename = "PACT_simulated_data.txt";
sensor_radius = 15e-3; % UST radius in [m]
num_of_A_lines = 100; % number of A-lines
sensor.frequency_response = [2.25e6 70]; % central frequency and bandwidth 
                                         % of the UST
sampling_frequency = 25;     % sampling frequency in [MHz]
number_of_time_points = 500; % number of time points in each A-line
object_sim.Nx = 400;  % number of grid points in the x (row) direction
object_sim.Ny = 400;  % number of grid points in the y (column) direction
object_sim.x = 40e-3; % total grid size [m]
object_sim.y = 40e-3; % total grid size [m]
%%%
time.dt = 1/(sampling_frequency*1e6); % sampling time in seconds
time.length = number_of_time_points;  % number of time points
time.t_array = 0:1:time.length-1;     % time array of time steps
time.t_array = time.t_array*time.dt;
Nx = object_sim.Nx;
Ny = object_sim.Ny;
dx = object_sim.x/object_sim.Nx; % grid point spacing in the x direction
dy = object_sim.y/object_sim.Ny; % grid point spacing in the y direction

kgrid = kWaveGrid(object_sim.Nx, dx, object_sim.Ny, dy);
kgrid.t_array = (0:time.length-1)*time.dt;
cart_sensor_mask = makeCartCircle(sensor_radius, num_of_A_lines);
sensor.mask = cart_sensor_mask;

%% Creating numerical point targets (3 points placed along the x axis)
source.p0 = zeros(object_sim.Nx, object_sim.Ny);
source.p0(205, 200) = 1; source.p0(205, 280) = 1; source.p0(205, 120) = 1;

%% set the input arguments for k wave
input_args = {"Smooth", false, "PMLInside", false, "PlotPML", ...
    false,'PlotSim', false};
%% Simulate the k-wave
simulated_PAT_data = kspaceFirstOrder2D(kgrid,medium,source,sensor, ...
    input_args{:});
%% plot the simulated A-line
figure, plot(simulated_PAT_data(50,:),'b',"LineWidth",2);
set(gca,"LineWidth",2,"XTick",[50 150 250 350 450],"YTick", ...
    [-4e-3 0 4e-3],"fontweight", "bold","fontsize",14);
title("A-line","fontweight","bold","fontsize",16);
xlabel("Time points","fontweight","bold", "fontsize",14);
ylabel("PA Signal amplitude","fontweight","bold", "fontsize",14);

%% plot PAT Sinogram
figure, imagesc(simulated_PAT_data);
set(gca,"LineWidth",2,"XTick",[50 150 250 350 450],"YTick", ...
    [10 30 50 70 90],"fontweight", "bold","fontsize",14);
title("PACT Sinogram","fontweigh","bold","fontsize",16);
xlabel("Time points","fontweight","bold", "fontsize",14);
ylabel("UST positions","fontweight","bold", "fontsize",14);

save(outputfilename, "simulated_PAT_data", "-ascii");
