clear all;
close all;

% Parameters
Ro = 1000;  % kg/m^3
C = 1500;   % m/s
freq = 20e6; 
lambda = C / freq; 
ka = 2 * pi / lambda;
R_foc = 9e-3; 
X0 = 0; Z0 = 0; Y0 = 0;
aperture_size = lambda * 2;
a = aperture_size;
Nsi = 200; 
N_phi_max = 200;
si_max = pi / 3;
element_Number = 63;
element_alfa = linspace(0, 2 * pi, element_Number + 1); 
element_fi = ones(1, element_Number + 1) * pi / 8;
element_ui = ones(1, element_Number + 1); 

% Modify element_ui to create a "T" shape
% Adjust elements' amplitudes to shape the letter "T"
% Central column of the "T"
element_ui(20:23) = 10;  % Increase amplitude for central elements
% Top bar of the "T"
element_ui(1:5) = 5;     % Increase amplitude for top elements
element_ui(59:63) = 5;   % Increase amplitude for top elements

Nsr = 10;

% Visualization of the transducer elements
figure; hold on;
[X, Y, Z] = sphere;
surf(0.8 * R_foc * X, 0.8 * R_foc * Y, 0.8 * R_foc * Z);

% Calculate pressure field
for i = 1:element_Number
    fi = element_fi(i);
    alfa = element_alfa(i);
    for y2 = -a/2:a/Nsr:a/2
        for z2 = -a/2:a/Nsr:a/2
            XYZ1 = [cos(fi) 0 sin(fi); 0 1 0; -sin(fi) 0 cos(fi)] * [-R_foc y2 z2]';
            XYZ = [1 0 0; 0 cos(alfa) -sin(alfa); 0 sin(alfa) cos(alfa)] * XYZ1;
            plot3(XYZ(1), XYZ(2), XYZ(3), '.b');
        end
    end
end

% Calculate intensity at different field points
i = 0;
for si = 0:(si_max) / (Nsi - 1):si_max
    i = i + 1;
    j = 0;
    Nphi(i) = max(round(N_phi_max * si / si_max + 1), 1);
    for phi = 0:(2 * pi - 0) / (Nphi(i) - 1):2 * pi
        j = j + 1;
        Xp(i, j) = R_foc * cos(si);
        Yp(i, j) = R_foc * sin(si) * cos(phi);
        Zp(i, j) = R_foc * sin(si) * sin(phi);
        fi_data(i, j) = phi;
        si_data(i, j) = si;
    end
end

% Calculate pressure and intensity
for i = 1:Nsi
    for j = 1:Nphi(i)
        P(i, j, 1:element_Number) = 0;
        for k = 1:element_Number
            fi = element_fi(k);
            alfa = element_alfa(k);
            XYZ1 = [1 0 0; 0 cos(alfa) sin(alfa); 0 -sin(alfa) cos(alfa)] * [Xp(i,j) Yp(i,j) Zp(i,j)]';
            XYZ2 = [cos(fi) 0 -sin(fi); 0 1 0; sin(fi) 0 cos(fi)] * XYZ1;
            XYZbti = XYZ2 + [R_foc 0 0]';
            P(i, j, k) = single_element_US_reduced(0, 0, aperture_size, element_ui(k), XYZbti(2), XYZbti(3), XYZbti(1));
        end
        I_p1(i, j) = (abs(sum(P(i, j, :))))^2;
    end
end

% Normalize intensity
I_p1 = I_p1 / max(max(I_p1));

% Plotting the results
figure; hold on;
Ipmin = min(min(I_p1));
Ipmax = max(max(I_p1));
for i = 1:Nsi
    for j = 1:Nphi(i)
        x = (I_p1(i,j) - Ipmin) / (Ipmax - Ipmin);
        x = min(max(x, 0), 1);
        plot3(Xp(i,j), Yp(i,j), Zp(i,j), '.', 'Color', [x 0 1-x], 'MarkerSize', 4);
    end
end

% Additional plotting for the T shape cross-section
fi_sec = [0 pi/3 pi/2 pi 4*pi/3 3*pi/2 2*pi];
figure; hold on; 
for j = 1:length(fi_sec)
    for i = 1:Nsi
        [val, ind] = min((fi_data(i, 1:Nphi(i)) - fi_sec(j)).^2);
        P_sec(i, j) = I_p1(i, ind);
    end
end
plot(P_sec);
