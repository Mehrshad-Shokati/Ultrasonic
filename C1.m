clear all;
close all;

disp('Starting the simulation...');

% Checkpoint 1: Initialization
disp('Checkpoint 1: Initialization...');
% Physical Constants
Ro = 1000;  % kg/m^3 (Density of the medium, e.g., water)
C = 1500;   % m/s (Speed of sound in the medium)
freq = 20e6; % 20 MHz (Frequency of the ultrasound wave)
lambda = C/freq; % Wavelength of the sound in meters
ka = 2*pi/lambda; % Wave number

disp('Constants initialized.');

% Checkpoint 2: Simulation Parameters
disp('Checkpoint 2: Setting Simulation Parameters...');
% Simulation Parameters
R_foc = 9e-3;  % Focal length in meters
aperture_size = lambda * 2; % Size of the aperture in meters
Nsi = 50; % Number of sampling points for 'si' (reduced to speed up testing)
N_phi_max = 50; % Maximum number of angular samples (reduced to speed up testing)
si_max = pi/3; % Maximum angle in radians
element_Number = 63; % Number of elements in the annular array

% Checkpoint 3: Array Setup
disp('Checkpoint 3: Array Setup...');
% Transducer Array Setup
element_alfa = linspace(0, 2*pi, element_Number + 1); % Angular positions for each array element
element_alfa = element_alfa(1:end-1); % Removing duplicate point at 2*pi
element_fi = pi/8 * ones(1, element_Number); % Tilt angle for each element
element_ui = zeros(1, element_Number); % Initialize all excitations to zero

% Define the indices to form the "C" pattern
C_indices = [5, 6, 7, 8, 9, 15, 22, 23, 24, 25, 26, 27, 28, 34, 40]; 
element_ui(C_indices) = 1.0;  % Strong excitation for the "C" shape

% Visualization Setup
Nsr = 10; % Number of subdivisions for element integration
[X, Y, Z] = sphere; % Creating a unit sphere for visualization
figure; hold on;
surf(0.8 * R_foc * X, 0.8 * R_foc * Y, 0.8 * R_foc * Z); % Displaying the focal sphere

% Plotting the array elements
for i = 1:element_Number
    fi = element_fi(i);
    alfa = element_alfa(i);
    for y2 = -aperture_size/2:aperture_size/Nsr:aperture_size/2
        for z2 = -aperture_size/2:aperture_size/Nsr:aperture_size/2
            XYZ1 = [cos(fi) 0 sin(fi); 0 1 0; -sin(fi) 0 cos(fi)] * [-R_foc y2 z2]';
            XYZ = [1 0 0; 0 cos(alfa) -sin(alfa); 0 sin(alfa) cos(alfa)] * XYZ1;
            plot3(XYZ(1), XYZ(2), XYZ(3), '.b');
        end
    end
end
disp('Checkpoint 3 Complete: Array setup and visualization done.');

% Checkpoint 4: Generating Pressure Field Points
disp('Checkpoint 4: Generating Pressure Field Points...');
% Generate the pressure field points on a spherical surface
Xp = zeros(Nsi, N_phi_max); 
Yp = zeros(Nsi, N_phi_max); 
Zp = zeros(Nsi, N_phi_max); 
I_p1 = zeros(Nsi, N_phi_max);

for i = 1:Nsi
    if mod(i, 10) == 0
        disp(['Generating pressure points, si index: ', num2str(i), '/', num2str(Nsi)]);
    end
    si = si_max * (i - 1) / (Nsi - 1);
    Nphi = max(round(N_phi_max * si / si_max + 1), 1);
    for j = 1:Nphi
        phi = 2 * pi * (j - 1) / (Nphi - 1);
        Xp(i, j) = R_foc * cos(si);
        Yp(i, j) = R_foc * sin(si) * cos(phi);
        Zp(i, j) = R_foc * sin(si) * sin(phi);
    end
end
disp('Checkpoint 4 Complete: Pressure field points generated.');

% Checkpoint 5: Calculating Pressure Field
disp('Checkpoint 5: Calculating Pressure Field...');
% Calculate the pressure field for each point on the spherical surface
counter = 0; % Initialize a counter
for i = 1:Nsi
    if mod(i, 10) == 0
        disp(['Calculating pressure field, si index: ', num2str(i), '/', num2str(Nsi)]);
    end
    for j = 1:N_phi_max
        P = zeros(1, element_Number);
        for k = 1:element_Number
            fi = element_fi(k);
            alfa = element_alfa(k);
            XYZ1 = [1 0 0; 0 cos(alfa) sin(alfa); 0 -sin(alfa) cos(alfa)] * [Xp(i,j) Yp(i,j) Zp(i,j)]';
            XYZ2 = [cos(fi) 0 -sin(fi); 0 1 0; sin(fi) 0 cos(fi)] * XYZ1;
            XYZbti = XYZ2 + [R_foc 0 0]';
            P(k) = single_element_US_reduced(0, 0, aperture_size, element_ui(k), XYZbti(2), XYZbti(3), XYZbti(1));
        end
        I_p1(i,j) = (abs(sum(P)))^2;
    end
end
disp('Checkpoint 5 Complete: Pressure field calculation done.');

% Checkpoint 6: Normalizing the Intensity Field
disp('Checkpoint 6: Normalizing the Intensity Field...');
% Normalize the intensity field
I_p1 = I_p1 / max(I_p1(:));
disp('Checkpoint 6 Complete: Intensity field normalized.');

% Ensure all matrices are the same size
minRows = min([size(Xp, 1), size(Yp, 1), size(I_p1, 1)]);
minCols = min([size(Xp, 2), size(Yp, 2), size(I_p1, 2)]);

Xp = Xp(1:minRows, 1:minCols);
Yp = Yp(1:minRows, 1:minCols);
I_p1 = I_p1(1:minRows, 1:minCols);

% Checkpoint 7: 2D Contour Visualization
disp('Checkpoint 7: 2D Contour Visualization...');
% 2D Contour Visualization with Depth
figure;
contourf(Xp, Yp, I_p1, 20, 'LineColor', 'none'); % 20 contour levels
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('Pressure Distribution in the X-Y Plane');
axis equal;
disp('Checkpoint 7 Complete: 2D contour plot done.');

% Function to calculate the pressure at a point due to a single element
function [P0] = single_element_US_reduced(X0, Y0, aperture_Lxy, un, Xp, Yp, Zp)
    % X0, Y0: Center of the source element
    % aperture_Lxy: Element length
    % Xp, Yp, Zp: Coordinates of the pressure field point
    
    % Constants and Parameters
    Ro = 1000;  % kg/m^3 (Density)
    C = 1500;   % m/s (Speed of sound)
    freq = 10e6; % 10 MHz
    lambda = C / freq; % Wavelength
    ka = 2 * pi / lambda; % Wave number
    alfa = 0; % Attenuation coefficient
    
    % Mesh and Integration Setup
    Nxy = 201; % Reduced mesh size for faster execution
    x0 = X0 + linspace(-0.5, 0.5, Nxy) * aperture_Lxy;
    y0 = Y0 + linspace(-0.5, 0.5, Nxy) * aperture_Lxy;
    h_xy = x0(2) - x0(1); % Grid spacing
    [X0_2D, Y0_2D] = meshgrid(x0, y0);
    
    % Calculate distance to the field point
    rho_p1 = sqrt((X0_2D - Xp).^2 + (Y0_2D - Yp).^2 + Zp^2);
    
    % Compute the pressure field contribution using Rayleigh-Sommerfeld method
    f_xy_RS = exp(-(alfa + ka * 1i) * rho_p1) ./ rho_p1;
    integral_2D_RS = integral2D(h_xy, Nxy, f_xy_RS);
    P0 = 1i * Ro * C / lambda * un * integral_2D_RS;
end

% Function to compute a 2D integral using the trapezoidal rule
function [integral_2D] = integral2D(h, N, f_xy)
    integral_2D = h^2 / 4 * (f_xy(1,1) + f_xy(1,N) + f_xy(N,1) + f_xy(N,N) + ...
        4 * sum(f_xy(:)) + ...
        2 * (sum(f_xy(1,:)) + sum(f_xy(N,:)) + sum(f_xy(:,1)) + sum(f_xy(:,N)) - ...
        2 * f_xy(1,1) - f_xy(1,N) - f_xy(N,1)));
end


