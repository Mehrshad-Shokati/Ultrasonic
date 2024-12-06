clear all;
close all;

disp('Starting the simulation...');

% Physical Constants
Ro = 1000;  % kg/m^3 (Density of the medium, e.g., water)
C = 1500;   % m/s (Speed of sound in the medium)
freq = 20e6; % 20 MHz (Frequency of the ultrasound wave)
lambda = C/freq; % Wavelength of the sound
ka = 2*pi/lambda; % Wave number

disp('Constants initialized.');

% Simulation Parameters
R_foc = 9e-3;  % Focal length (distance from the array to the focal point)
aperture_size = lambda * 2; % Size of the aperture
Nsi = 200; % Number of sampling points for 'si'
N_phi_max = 200; % Maximum number of angular samples
si_max = pi/3; % Maximum angle for the simulation
element_Number = 63; % Number of elements in the annular array

disp('Simulation parameters set.');

% Transducer Array Setup
element_alfa = linspace(0, 2*pi, element_Number + 1); % Angular positions for each array element
element_alfa = element_alfa(1:end-1); % Remove duplicate point at 2*pi
element_fi = pi/8 * ones(1, element_Number); % Tilt angle for each element

% Modify the excitation for each element to form a "T" pattern
element_ui = zeros(1, element_Number); % Initialize all excitations to zero

% Define the indices for the vertical part of the "T"
vertical_indices = [29, 30, 31, 32, 33]; 

% Define the indices for the horizontal part of the "T"
horizontal_indices = [16, 30, 46]; 

% Set the excitation for the "T" pattern
element_ui(vertical_indices) = 1.0;  % Strong excitation for the vertical part
element_ui(horizontal_indices) = 1.0; % Strong excitation for the horizontal part

disp('Element excitation set.');

% Number of subdivisions for element integration
Nsr = 10;

% Initialize Pressure Field and Points for Visualization
Xp = zeros(Nsi, N_phi_max); 
Yp = zeros(Nsi, N_phi_max); 
Zp = zeros(Nsi, N_phi_max); 
I_p1 = zeros(Nsi, N_phi_max);

disp('Pressure field initialization complete.');

% Generate the pressure field points on a spherical surface
for i = 1:Nsi
    si = si_max * (i - 1) / (Nsi - 1);
    Nphi = max(round(N_phi_max * si / si_max + 1), 1);
    for j = 1:Nphi
        phi = 2 * pi * (j - 1) / (Nphi - 1);
        Xp(i,j) = R_foc * cos(si);
        Yp(i,j) = R_foc * sin(si) * cos(phi);
        Zp(i,j) = R_foc * sin(si) * sin(phi);
    end
end

disp('Pressure field points generated.');

% Calculate the pressure field for each point on the spherical surface
for i = 1:Nsi
    for j = 1:Nphi
        P = zeros(1, element_Number);
        for k = 1:element_Number
            fi = element_fi(k);
            alfa = element_alfa(k);
            XYZ1 = [1 0 0; 0 cos(alfa) sin(alfa); 0 -sin(alfa) cos(alfa)] * [Xp(i,j) Yp(i,j) Zp(i,j)]';
            XYZ2 = [cos(fi) 0 -sin(fi); 0 1 0; sin(fi) 0 cos(fi)] * XYZ1;
            XYZbti = XYZ2 + [R_foc 0 0]';
            P(k) = single_element_US_ReySumerfield(0, 0, aperture_size, element_ui(k), XYZbti(2), XYZbti(3), XYZbti(1));
        end
        I_p1(i,j) = (abs(sum(P)))^2;
    end
end

disp('Pressure field calculation complete.');

% Normalize the intensity field
I_p1 = I_p1 / max(I_p1(:));

disp('Intensity field normalized.');

% Basic Check for a Simple Plot
figure;
plot(Xp(:), Yp(:), '.');
title('Simple XY Plot');
xlabel('X');
ylabel('Y');
grid on;

disp('Simple XY plot done.');

% Visualization of the Pressure Field (3D Surface Plot)
figure;
surf(Xp, Yp, I_p1, 'EdgeColor', 'none');
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Acoustic Pressure (Pa)');
title('3D Visualization of the Ultrasound Field');
view(2); % To view from top
axis tight;

disp('3D surface plot done.');

% 2D Contour Visualization
figure;
contourf(Xp, Yp, I_p1, 20, 'LineColor', 'none'); % 20 contour levels
colorbar;
xlabel('X (mm)');
ylabel('Y (mm)');
title('Pressure Distribution in the X-Y Plane');
axis equal;

disp('2D contour plot done.');

% Function to calculate the pressure at a point due to a single element
function [P0] = single_element_US_ReySumerfield(X0, Y0, aperture_Lxy, un, Xp, Yp, Zp)
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
    Nxy = 401;
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
