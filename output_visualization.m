output = readmatrix('output\FOUS_output_dt=0.0100_dx=0.0300.csv');   % 根据需要修改文件名(FOUS/LWS/UWBS和dt、dx)

[T, N] = size(output);
x = linspace(0, 3, N);
dt = 0.01;

figure;
hold on;

% 计算并绘制初始时刻的参考函数（红色）
t_current = 1*dt;
u_ref = sin(2*pi*(x - t_current));
h_ref = plot(x, u_ref, 'r-');

h_sim = plot(x, output(1, :), 'b-');  % 数值解图形句柄h_sim（蓝色）

xlabel('x');
ylabel('u');
title(sprintf('时间: %.4f', t_current));  %标题
legend('数值解', '解析解', 'Location', 'northeast');  %标注
axis([0 3 -1.2 1.2]);        % 固定坐标轴范围
grid on;

for t = 1:T
    set(h_sim, 'YData', output(t, :));  % 更新数值解
    
    t_current = t*dt;
    u_ref = sin(2*pi*(x - t_current));
    set(h_ref, 'YData', u_ref);    % 更新参考函数
    
    title(sprintf('时间: %.4f', t_current));  % 更新标题
    
    pause(0.01);  % 控制动画速度
    drawnow;
end
hold off;