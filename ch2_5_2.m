% This program is written to simulate the effect of an finite pulse width 
% Gaussian excitation (1 us pulse width) of a spherical object of radius
% 0.2 mm using K wave. The photoacoustic wave generated from the object is  
% captured by a point sensor (detector) kept at a distance of 2 mm from the  
% centre of the object. 

% author: Manojit Pramanik
% date: 05 August 2026
% k-Wave Version 1.4.1 and Matlab R2026a Update 4 (26.1.0.3312084)

clear all;

% =========================================================================
% SIMULATION
% =========================================================================
% defining the geometry of computations
sensor_radius = 2e-3; %[m]  
ball_radius = 0.2e-3; %[m]
pulsewidth_min = 0.1e-6; % min width of the pulse which excites the source.
pulse_width = 1e-6; % Gaussian pulse width

% create the computational grid
dx = 0.025e-3;    	% grid point spacing in the x direction [m]
dy = dx;            % grid point spacing in the y direction [m]
dz = dx;            % grid point spacing in the z direction [m]
PMLe = 25*dx; %[no of grid points for PML*dx]
grig_length = (2*sensor_radius)+(2*PMLe); % total grid length in a
                                          % particular direction [m]
Nx = grig_length/dx;    % number of grid points in the x (row) direction
Ny = grig_length/dx;    % number of grid points in the y (column) direction
Nz = grig_length/dx;    % number of grid points in the z (axis) direction

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

% define the properties of the propagation medium
medium.sound_speed = 1500;  % [m/s]

% create the time array
time.dt = pulsewidth_min/10;  % sampling interval[s]
time.total =  1.5e-6 + sensor_radius/medium.sound_speed; % Total simulation
                 % time or the time taken by the wave to reach the detector
time.Nt = ceil(time.total/time.dt); % No of time steps[points]
time.t_array = 0:1:time.Nt-1; % time array having Nt time steps
time.t_array = time.t_array*time.dt; % array containing the values of time.
kgrid.t_array = time.t_array;

% create initial pressure distribution using makeBall
ball_magnitude = 1;       % [au]
ball_x_pos = Nx/2;        % [grid points]
ball_y_pos = Ny/2;        % [grid points]
ball_z_pos = Nz/2;        % [grid points]
ball_radius = ball_radius/dx; % [grid points]
ball_1 = ball_magnitude*makeBall(Nx, Ny, Nz, ball_x_pos, ball_y_pos, ...
    ball_z_pos, ball_radius);
source.p_mask = ball_1;

% defining the excitation pulse type - Gaussian
Np = pulse_width/time.dt; % Number of time points of pulse signal
mu=Np/2;

source_gaussianpulse = zeros(1,time.Nt);
sigma = Np/6; %Variance of guassian pdf
%mean = 0;
for j=1:Np
        exponent =(((j-mu)^2)/((2*sigma^2)));
    %value =(1/(sqrt((2*pi*sigma^2))))* (exp(-exponent));
    value = (exp(-exponent));
    source_gaussianpulse(j) = value;
end
for j=Np+1:time.Nt
    source_gaussianpulse(j)=0;
end
source.p =source_gaussianpulse;
kgrid.t_array = time.t_array;

% smooth the source
% source.p = filterTimeSeries(kgrid, medium, source.p);

sensor_data = zeros(1,time.Nt); % To record the PA data

% define a single sensor element (point sensor)
sensor_distance = sensor_radius / dx ;
sensor.mask = zeros(Nx, Ny, Nz);
sensor.mask(Nx/2 +(sensor_distance), Ny/2, Nz/2) = 1;

% run the simulation in 3D
sensor_data(:)= kspaceFirstOrder3D(kgrid, medium, source, sensor);
% =========================================================================
% VISUALISATION
% =========================================================================
% plot the simulated sensor data
figure;
[t_sc, scale, prefix] = scaleSI(max(time.t_array(:)));

%plotting a portion of time
t1 = 1;
t2 = time.Nt;
t_array1 = time.t_array*scale;
t_plot = t_array1(t1:t2);
sensordata1=sensor_data(:);
sensorplot1=sensordata1(t1:t2);
% normalisedsensorplot1=(sensorplot1)./(max(sensorplot1));
normalisedsensorplot1=sensorplot1;
hold on

% plot for 0.2 mm source radius
plot(t_plot, normalisedsensorplot1, 'b-','LineWidth',2.5); %blue
tick = [0.01 0.01];
set(0, 'DefaultAxesTickLength', tick);

ax = gca;
ax.FontSize = 14;                    % tick labels and axis numbers
ax.FontWeight = 'normal';            % optional: 'bold'
ax.LineWidth = 1;                  % axis box / tick line thickness

xlabel(['Time [µs]'],'FontSize',14);
ylabel('Signal Amplitude','FontSize',14);
axis tight;

% plotting the Gaussian pulse
figure,
plot(t_plot, source_gaussianpulse, 'r-','LineWidth',2.5);
tick = [0.01 0.01];
set(0, 'DefaultAxesTickLength', tick);

ax = gca;
ax.FontSize = 14;                    % tick labels and axis numbers
ax.FontWeight = 'normal';            % optional: 'bold'
ax.LineWidth = 1;                  % axis box / tick line thickness

xlabel(['Time [µs]'],'FontSize',14);
ylabel('Signal Amplitude','FontSize',14);
axis tight;

save ch2_5_2.mat
