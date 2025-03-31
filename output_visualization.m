output = readmatrix('output\FOUS_output_dt=0.0100_dx=0.0300.csv');   % ������Ҫ�޸��ļ���(FOUS/LWS/UWBS��dt��dx)

[T, N] = size(output);
x = linspace(0, 3, N);
dt = 0.01;

figure;
hold on;

% ���㲢���Ƴ�ʼʱ�̵Ĳο���������ɫ��
t_current = 1*dt;
u_ref = sin(2*pi*(x - t_current));
h_ref = plot(x, u_ref, 'r-');

h_sim = plot(x, output(1, :), 'b-');  % ��ֵ��ͼ�ξ��h_sim����ɫ��

xlabel('x');
ylabel('u');
title(sprintf('ʱ��: %.4f', t_current));  %����
legend('��ֵ��', '������', 'Location', 'northeast');  %��ע
axis([0 3 -1.2 1.2]);        % �̶������᷶Χ
grid on;

for t = 1:T
    set(h_sim, 'YData', output(t, :));  % ������ֵ��
    
    t_current = t*dt;
    u_ref = sin(2*pi*(x - t_current));
    set(h_ref, 'YData', u_ref);    % ���²ο�����
    
    title(sprintf('ʱ��: %.4f', t_current));  % ���±���
    
    pause(0.01);  % ���ƶ����ٶ�
    drawnow;
end
hold off;