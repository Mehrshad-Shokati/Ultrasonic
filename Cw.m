% ایجاد شبکه برای شبکیه چشم  
theta = linspace(0, 2*pi, 100); % زاویه  
r = linspace(0, 5, 50); % شعاع  

% ایجاد مختصات قطبی برای شبکیه  
[R, T] = meshgrid(r, theta);  
[X, Y] = pol2cart(T, R);  

% مدل سه‌بعدی شبکیه چشم  
Z = sin(R) .* cos(T); % تابع دلخواه برای نمایش شکل شبکیه  

% تعریف ترانسدیوسرها  
num_transducers = 8; % تعداد ترانسدیوسرها  
transducer_angles = linspace(0, 2*pi, num_transducers+1);  
transducer_angles(end) = []; % حذف آخرین زاویه چون تکراری است  
transducer_radii = 3; % شعاع ترانسدیوسرها  

% اعمال الگوی C  
pattern_C = zeros(size(Z));  
for i = 1:num_transducers  
    % مختصات ترانسدیوسر محاسبه می‌شود  
    x_t = transducer_radii * cos(transducer_angles(i));  
    y_t = transducer_radii * sin(transducer_angles(i));  
    
    % ایجاد یک الگو حول هر ترانسدیوسر (مدل ساده)  
    distance_from_transducer = sqrt((X - x_t).^2 + (Y - y_t).^2);  
    pattern_C = pattern_C + exp(-distance_from_transducer.^2); % تابع گوسی  
end  

% ترکیب الگوی C با شبکیه  
Z_combined = Z + pattern_C;   

% ترسیم مدل شبکیه همراه با الگوی C  
figure;  
surf(X, Y, Z_combined, 'EdgeColor', 'none'); % استفاده از surf برای ترسیم  
colormap(jet); % تغییر رنگ به تابع رنگی جت  
colorbar; % نمایش نوار رنگ  
xlabel('X coordinate');  
ylabel('Y coordinate');  
zlabel('Z coordinate');  
title('3D Model of Retina with Pattern C');  
axis equal; % تنظیم مقیاس محور‌ها برابر  
view(30, 30); % زاویه دید مناسب  
grid on; % نمایش شبکه