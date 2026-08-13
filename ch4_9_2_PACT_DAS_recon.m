% PACT Delay-and-sum reconstruction for circular scan with simulated data
% generated earlier

% author: Manojit Pramanik
% date: 05 August 2026
% k-Wave Version 1.4.1 and Matlab R2026a Update 4 (26.1.0.3312084)

clc; clear all; close all;
%% —————— parameters ————————————
filename_pa = "PACT_simulated_data.txt"; % input data
fs = 25; % sampling frequency (in MHz)
x_size = 30e-3; y_size = 30e-3; % reconstructed image size in m
Npixel_x = 300; Npixel_y = 300; % number of  pixels in x and y directions
radius = 250;  % reconstruction radius
soundv = 1500; % velocity of sound in m/s

pa_data = load(filename_pa).';        % loading the input data
[Nsample, Nposition] = size(pa_data); % Nposition: number of UST positions
                                      % Nsample: number of time points
angle_per_step = 2*pi/Nposition;
R = radius/(fs*1e6)*soundv; % conversion of reconstruction radius in to m

%% generating the reconstruction grid
x_scale = Npixel_x/x_size; y_scale = Npixel_y/y_size;
x_img = ((1:Npixel_x)-((Npixel_x+1)/2))*x_size/(Npixel_x-1);
x_img = ones(Npixel_y,1)*x_img;
y_img = (((Npixel_y+1)/2)-(1:Npixel_y))*y_size/(Npixel_y-1);
y_img = y_img.'*ones(1,Npixel_x);
disp("Please keep patience, PACT image reconstruction is going on. . . ");

%% Delay-and-sum
x_receive = cos((0:Nposition-1)*angle_per_step)*R; % UST x coordinates
y_receive = sin((0:Nposition-1)*angle_per_step)*R; % UST y coordinates
pa_img = zeros(Npixel_x, Npixel_y);

for i = 1:Nposition
    pa_data_tmp = [pa_data(:,i); 0];
    dist = sqrt((x_img-x_receive(i)).^2+(y_img-y_receive(i)).^2);
    % distance of the each pixel in the grid for each UST position
    pxd = round(dist/soundv*fs*1e6);
    % index
    inrange = (pxd>=1)&(pxd<=Nsample)&(sqrt(x_img.*x_img+y_img.*y_img)< R);
    pxd = (inrange).*pxd+(1-inrange).*(Nsample+1);
    pa_img = pa_img+pa_data_tmp(pxd);
end

%% Displaying the reconstructed image
X = linspace(1,Npixel_x,Npixel_x)*x_size/Npixel_x;
Y = linspace(1,Npixel_y,Npixel_y)*y_size/Npixel_y;
figure; imagesc(X*1e3,Y*1e3,fliplr(pa_img)); axis image;
value = [1:-1/63:0].'; gray2 = [value value value]; colormap(gray2);
set(gca,"LineWidth",1,"XTick",[5 15 25 35],"YTick",[5 15 25 35], ...
    "fontweight","bold","fontsize",12);
xlabel("mm","fontweight","bold", "fontsize",14);
ylabel("mm","fontweight","bold", "fontsize",15);
title("Reconstructed PACT image","fontweight","bold", "fontsize",16);
