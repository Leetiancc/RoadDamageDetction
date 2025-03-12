% 示例输入数据
clear;
clc;
eta = [15,15.5,15.4,15.3,15.2,15.1,14.5,14.6,14.7,14.8,14.9,13.5,13.6,13.7,13.8,13.9];
% eta = [12,12,12,12,11,12,11,12,12,11];
x0 = [0.5, 15]; % 初始猜测值 [γ, δ]

% 调用 fsolve
options = optimoptions('fsolve', 'Display', 'iter', 'TolFun', 1e-12, 'TolX', 1e-12);
[x, fval, exitflag] = fsolve(@(x) myNonlinearSystem(x, eta), x0, options);


% 输出结果
gamma = x(1);
delta = x(2);
disp(['gamma = ', num2str(gamma)]);
disp(['delta = ', num2str(delta)]);

x = 0:0.1:200; % x 从 0 到 20

% 计算 Gamma 分布的概率密度
y = gampdf(x, delta, gamma);
figure;
plot(x, y, 'b-', 'LineWidth', 2);
xlabel('x');
ylabel('Probability Density');
title(['Gamma Distribution with k=', num2str(delta), ', \theta=', num2str(gamma)]);
grid on;

%%
clear all;
clc;
for i = 1:15
    eta(i) = 30*rand;
end
x1 = 0.1:0.1:10;
x2 = 1:1:100;
for i = 1:length(x1)
    for j = 1:length(x2)
        gamma = x1(i);
        delta = x2(j);
        F = myNonlinearSystem([gamma, delta], eta);
        y1(i,j) = F(1);
        y2(i,j) = F(2);
    end
end

[X1, X2] = meshgrid(x1, x2);

figure
surf(x1,x2,y1,'FaceColor','interp','EdgeColor','none');
colorbar;
hold on;

contour_data1 = contourc(x1, x2, y1, [0,0]); % 自定义等高线
% 提取等高线坐标
x_contour = contour_data1(1, 2:end); % x 坐标
y_contour = contour_data1(2, 2:end); % y 坐标

contour_data2 = contourc(x1, x2, y2, [0,0]); % 自定义等高线
% 提取等高线坐标
x_contour2 = contour_data2(1, 2:end); % x 坐标
y_contour2 = contour_data2(2, 2:end); % y 坐标

% 绘制提取的等高线
figure;
plot(x_contour, y_contour, 'b-', 'LineWidth', 1); % 使用 plot 绘制
hold on
plot(x_contour2, y_contour2, 'b-', 'LineWidth', 1); % 使用 plot 绘制
xlabel('x');
ylabel('y');
title('Extracted Contour for z = 0');
grid on;
%%
figure;
surf(x1,x2,y2,'FaceColor','interp','EdgeColor','none');
colorbar;
hold on;
contour3(x1, x2, y2, 0, 'r', 'LineWidth', 1.5); % 自定义等高线





%% 函数
function F = myNonlinearSystem(x, eta)
gamma = x(1); % 第一个变量 γ
delta = x(2); % 第二个变量 δ
n = length(eta);

% 第一方程
F(1) = sum(log(gamma) + log(eta) - (psi(delta)));
% 第二方程
F(2) = sum(delta / gamma - eta);
end