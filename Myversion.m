% Main Script for Ultrasonic Field Simulation
% This code simulates the pressure field on the retina due to a US annular array.
% The stimulation of each transducer element can be adjusted to produce different output patterns.

clear all;
close all;

% Physical Constants
Ro = 1000;  % kg/m^3 (Density of the medium, e.g., water)
C = 1500;   % m/s (Speed of sound in the medium)
freq = 20e6; % 20 MHz (Frequency of the ultrasound wave)
lambda = C/freq; % Wavelength of the sound
ka = 2*pi/lambda; % Wave number

% Simulation Parameters
R_foc = 9e-3;  % Focal length (distance from the array to the focal point)
aperture_size = lambda * 2; % Size of the aperture
Nsi = 200; % Number of sampling points for 'si'
N_phi_max = 200; % Maximum number of angular samples
si_max = pi/3; % Maximum angle for the simulation
element_Number = 63; % Number of elements in the annular array

% Transducer Array Setup
element_alfa = linspace(0, 2*pi, element_Number + 1); % Angular positions for each array element
element_alfa = element_alfa(1:end-1); % Remove duplicate point at 2*pi
element_fi = pi/8 * ones(1, element_Number); % Tilt angle for each element
element_ui = ones(1, element_Number); % Excitation amplitude for each element (modifiable)

% Number of subdivisions for element integration
Nsr = 10;

% Initialize Pressure Field and Points for Visualization
Xp = zeros(Nsi, N_phi_max); 
Yp = zeros(Nsi, N_phi_max); 
Zp = zeros(Nsi, N_phi_max); 
I_p1 = zeros(Nsi, N_phi_max);

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

% Normalize the intensity field
I_p1 = I_p1 / max(I_p1(:));

% Visualization of the Pressure Field (2D Contour Plot)
figure; hold on;
for i = 1:Nsi
    for j = 1:Nphi
        x = (I_p1(i,j) - min(I_p1(:))) / (max(I_p1(:)) - min(I_p1(:)));
        plot3(Xp(i,j), Yp(i,j), Zp(i,j), '.', 'Color', [x 0 1-x], 'MarkerSize', 4);
    end
end
xlabel('X (mm)');
ylabel('Y (mm)');
zlabel('Z (mm)');
title('3D Visualization of the Pressure Field');

% % 2D Contour Visualization
% figure;
% contourf(Yp(1:Nsi,1:Nphi), Zp(1:Nsi,1:Nphi), I_p1(1:Nsi,1:Nphi), 'LineColor', 'none');
% colorbar;
% xlabel('Y (mm)');
% ylabel('Z (mm)');
% title('2D Contour Plot of the Pressure Field');
% axis equal;
% 
% % Function to calculate the pressure at a point due to a single element
% function [P0] = single_element_US_ReySumerfield(X0, Y0, aperture_Lxy, un, Xp, Yp, Zp)
%     % X0, Y0: Center of the source element
%     % aperture_Lxy: Element length
%     % Xp, Yp, Zp: Coordinates of the pressure field point
%     
%     % Constants and Parameters
%     Ro = 1000;  % kg/m^3 (Density)
%     C = 1500;   % m/s (Speed of sound)
%     freq = 10e6; % 10 MHz
%     lambda = C / freq; % Wavelength
%     ka = 2 * pi / lambda; % Wave number
%     alfa = 0; % Attenuation coefficient
%     
%     % Mesh and Integration Setup
%     Nxy = 401;
%     x0 = X0 + linspace(-0.5, 0.5, Nxy) * aperture_Lxy;
%     y0 = Y0 + linspace(-0.5, 0.5, Nxy) * aperture_Lxy;
%     h_xy = x0(2) - x0(1); % Grid spacing
%     [X0_2D, Y0_2D] = meshgrid(x0, y0);
%     
%     % Calculate distance to the field point
%     rho_p1 = sqrt((X0_2D - Xp).^2 + (Y0_2D - Yp).^2 + Zp^2);
%     
%     % Compute the pressure field contribution using Rayleigh-Sommerfeld method
%     f_xy_RS = exp(-(alfa + ka * 1i) * rho_p1) ./ rho_p1;
%     integral_2D_RS = integral2D(h_xy, Nxy, f_xy_RS);
%     P0 = 1i * Ro * C / lambda * un * integral_2D_RS;
% end
% 
% % Function to compute a 2D integral using the trapezoidal rule
% function [integral_2D] = integral2D(h, N, f_xy)
%     integral_2D = h^2 / 4 * (f_xy(1,1) + f_xy(1,N) + f_xy(N,1) + f_xy(N,N) + ...
%         4 * sum(f_xy(:)) + ...
%         2 * (sum(f_xy(1,:)) + sum(f_xy(N,:)) + sum(f_xy(:,1)) + sum(f_xy(:,N)) - ...
%         2 * f_xy(1,1) - f_xy(1,N) - f_xy(N,1)));
% end
