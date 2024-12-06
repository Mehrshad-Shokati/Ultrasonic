% تعریف پارامترهای ثابت  
Ro = 1000;         % چگالی آب (kg/m^3)  
C = 1500;         % سرعت صوت در آب (m/s)  
freq = 50000;     % فرکانس (Hz)  
lambda = C / freq; % طول موج  
ka = 2 * pi / lambda; % عدد موج  

% تنظیمات منابع  
R_foc = 50;       % فاصله کانونی (m)  
num_elements = 20; % تعداد عناصر در آرایه  

% تنظیم آرایه مبدل  
angles = linspace(-pi/4, pi/4, num_elements); % زاویه عناصر  
positions = [R_foc * cos(angles)', R_foc * sin(angles)']; % موقعیت عناصر  

% تولید مصور از حرف "C"  
figure;  
hold on;  

% ترسیم حرف "C" با نقاط  
theta = linspace(-pi/2, pi/2, 100); % زاویه برای نیم دایره  
x_C = 10 * cos(theta); % مختصات X  
y_C = 10 * sin(theta); % مختصات Y  

plot(x_C, y_C, 'b', 'LineWidth', 2); % رسم رنگ آبی برای "C"  

% ترسیم موقعیت اعضاء  
plot(positions(:, 1), positions(:, 2), 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r'); % رسم دایره‌های قرمز برای مبدل‌ها  

axis equal;  
xlabel('X coordinate (m)');  
ylabel('Y coordinate (m)');  
title('Pattern Representation of Letter "C"');  
grid on;  
hold off;  

% پردازش و محاسبه فشار (در صورت لزوم)  
% این بخش بستگی به نیاز شما دارد که در اینجا نشان داده نشده است.