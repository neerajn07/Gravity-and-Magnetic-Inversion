% Forward Modeling of Gravity Data using Prism
% Author: Neeraj Nainwal
% Institution: Queen's University

% Clear workspace and command window, close all figures
clear all; clc;
close all;

% Gravitational constant in SI units (m^3 kg^-1 s^-2)
G = 6.67430 * 10^-11;

%% Define Parameters
cubeLength = 4;               % Length of the cube in meters
drho = 2000.0;                % Density contrast in kg/m^3
yp = 0.0;                     % Y-coordinate of the prism center
zp = [0 5 10];                % Z-coordinates of the prism center
xp = -20:1:20;                % X-coordinates of observation points

% Compute intermediate variables for the prism dimensions
dx = [-cubeLength/2, cubeLength/2];
dy = [-cubeLength/2, cubeLength/2];
dz = [-3, -3 - cubeLength/2];

% Preallocate matrices for results
rhog_z = zeros(8, length(xp));
summed_columns = zeros(1, length(xp));

%% Create Repetitive Cubes
numCubes = 1;          % Number of repetitive cubes
cubeSpacing = 1;       % Spacing between cubes in meters

for ii = 1:length(zp)
    for ix = 1:length(xp)
        count = 0;
        rhog_ix = zeros(8, 1);  
        for nc = 1:numCubes
            xShift = (nc - 1) * cubeSpacing;        
            for i = 1:2
                rhoxi = dx(i) - xp(ix) - xShift;      
                for j = 1:2
                    rhoyj = dy(j) - yp;     
                    for k = 1:2
                        count = count + 1;
                        rhozk = dz(k) - zp(ii);         
                        R_ijk = sqrt(rhozk^2 + rhoyj^2 + rhoxi^2);
                        theta = (rhoxi * rhoyj) / (rhozk * R_ijk);
                        muijk = (-1)^(i+j+k);             
                        % Calculate the gravitational effect
                        rhog_z(count, ix) = G * drho * muijk * ((rhozk * atan(theta)) ...
                            - (rhoxi * log(R_ijk + rhoyj)) ...
                            - (rhoyj * log(R_ijk + rhoxi)));         
                        rhog_ix(count) = rhog_z(count, ix);
                    end
                end
            end
        end
        summed_columns(ii, ix) = sum(rhog_ix);
    end
end

%% Data Plotting
% Plot the forward model results
subplot(2, 1, 1)
plot(xp, summed_columns(1,:), '-o', xp, summed_columns(2,:), '-o', xp, summed_columns(3,:), '-o', 'LineWidth', 1.5)
labels = num2str(zp(1:3).', 'Elevation %d m');
legend(labels, 'location', 'best') 
title('Forward Model')
ylabel('Gravity Anomaly (mGal)')
xlabel('Observation Points (km)')
hold on

% Plot the 3D model of the prism
subplot(2, 1, 2)
% Create cube vertices
vert = [-2 2 -3; 2 2 -3; -2 -2 -3; 2 -2 -3; -2 -2 -7; 2 2 -7; -2 2 -7; 2 -2 -7];
fac = [1 2 4 3; 2 6 8 4; 2 6 7 1; 1 7 5 3; 3 5 8 4; 5 8 6 7];

view(3)
axis equal
xlim([xp(1), xp(end)]);  
ylim([xp(1)/2, xp(end)/2]);  
zlim([xp(1), 0]);  
xlabel('X'); ylabel('Y'); zlabel('Z');
title('Model'); axis vis3d; grid on;

% Plot repetitive cubes
for nc = 1:numCubes
    xShift = (nc - 1) * cubeSpacing;  % Shift in x-coordinate for repetitive cubes
    vertShifted = vert;
    vertShifted(:, 1) = vertShifted(:, 1) + xShift;   
    patch('Vertices', vertShifted, 'Faces', fac, ...
        'FaceVertexCData', hsv(6), 'FaceColor', 'flat')
    hold on
end

% Save the plot
outputPath = fullfile('./results/', 'Forward_model.jpg');
saveas(gcf, outputPath);
