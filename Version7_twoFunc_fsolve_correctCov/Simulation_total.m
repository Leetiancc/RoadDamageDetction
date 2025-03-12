clear all;
clc;

tic;

Lambda_hat1 = [];
a1_lambdaSet_lo = [];
a1_lambdaSet_up = [];
a1_BootstrapSet = [];

data = readmatrix('GasLaser_5.csv');
m = length(data(:,2))./17;
Y = cell(m,1);
Fix_time = cell(m,1);
Cox_z = data(:,4);


% % col_min = min(Cox_z, [], 1)-1e-1;  % 1 x p
% % col_max = max(Cox_z, [], 1)+1e-1;  % 1 x p
% % 
% % 用 bsxfun 或者新版本 MATLAB 的自动广播
% % Cox_z = (Cox_z - col_min) ./ (col_max - col_min);
% % 
% % Cox_z = zscore(Cox_z);

for i = 1:m
    Y{i} = data(17*(i-1)+1:17*i,1);
    y = Y{i};
    y = y(2:end);
    Y{i} = y';
    Fix_time{i} = data(17*(i-1)+1:17*i,3);
    tt = Fix_time{i};
    tt = tt(2:end);
    Fix_time{i} = tt';
end
fix_time = 250:250:4000;

% figure;
% hold on;
% for i = 1:m
%     plot([0, Fix_time{i}], [0, Y{i}], '-o', 'MarkerSize', 4, 'MarkerFaceColor', 'b', 'LineWidth', 1);
% end
% yline(6, 'k-', 'linewidth', 1);
% xlim([0, 4000]);
% ylim([0, 12])
% title('GaslaserDegradation');
% xlabel('Time');
% ylabel('Electricity');

fprintf('模拟开始\n')

% 预设均值函数
psi1 = @(t) 0.002 * t; % 线性均值函数

D = 7; % 故障阈值

% 设定容忍度和最大迭代次数
tol = 1e-6;
max_iter = 1e3;
%%
% 计算lambda(t)系数初始值、eta初始值
[init_Lambda_hat1, init_lambda_hat1, init_delta, init_gamma, init_eta1, init_a1, init_beta] = initial_Lambda(Y, m, 'linear', fix_time, 10, 1 , 1e-4, Cox_z);
%
% 执行EM算法
[Lambda_hat1, lambda_hat1, a1, a1_eta, delta, gamma, beta, abcd] = EM_algorithm_hat1(Y, Fix_time, init_a1, tol, max_iter, 'linear', m, init_eta1, fix_time, init_lambda_hat1, init_Lambda_hat1, init_delta, init_gamma, init_beta, Cox_z);
std(a1_eta);
%%
% 
% t = linspace(0, 4000, 17); % 统一时间点以便于绘图
% %Bootstrap 方法计算 95% 置信区间
% %执行bootstrap函数
% 
% B = 10;
% [a1_BootstrapSet, a1_LambdahatSet, a1_lambdahatSet, eta_1bs, delta_1bs, gamma_1bs] = Bootstrap(m, fix_time, a1, a1_eta, 'linear', tol, max_iter, lambda_hat1, B, Fix_time, Y, Lambda_hat1, delta, gamma);
% 
% %计算置信区间：基于百分位数
% alpha = 0.025;
% a1_lo = quantile(a1_BootstrapSet, alpha);  % 2.5%分位数
% [~,idx] = min(abs(a1_BootstrapSet - a1_lo));
% a1_lambdaSet_lo = a1_LambdahatSet{idx};
% a1_up = quantile(a1_BootstrapSet, 1 - alpha);  % 97.5%分位数;
% [~,idx] = min(abs(a1_BootstrapSet - a1_up));
% a1_lambdaSet_up = a1_LambdahatSet{idx};
% delta_lo = quantile(delta_1bs, alpha);
% delta_up = quantile(delta_1bs, 1-alpha);
% gamma_lo = quantile(gamma_1bs, alpha);
% gamma_up = quantile(gamma_1bs, 1-alpha);
% 
% % eta_1lo = quantile(eta_1bs, alpha);  % 2.5%分位数
% % eta_1up = quantile(eta_1bs, 1 - alpha);  % 97.5%分位数;
% 
% % 绘制
% figure;
% subplot(1, 2, 1);
% hold on;
% for i = 1:m
%     plot([0,Fix_time{i}], [0,Y{i}], 'k'); % 绘制模拟退化路径
% end
% yline(D, 'k-', 'LineWidth', 1); % 故障阈值 D = 7
% xlabel('time');
% ylabel('degradation');
% title('Degradation paths (Linear)');
% 
% subplot(1, 2, 2);
% plot(t, [0,Lambda_hat1], 'k', 'LineWidth', 1.5); hold on;
% plot(t, psi1(t), 'k--', 'LineWidth', 1.5); % 真实均值函数
% plot(t, [0, a1_lambdaSet_lo], 'k:', 'LineWidth', 1.5); % 置信区间下限
% plot(t, [0, a1_lambdaSet_up], 'k:', 'LineWidth', 1.5); % 置信区间上限
% xlabel('time');
% ylabel(' \Lambda (t)');
% title('MLE and 95% CI (Linear)');
% 
% % 设置整体图形属性
% set(gcf, 'Position', [100, 100, 1500, 500]);


elapsed_time = toc;
disp(['代码运行时间为: ', num2str(elapsed_time), ' 秒']);

% pd = fitdist(init_eta1', 'Gamma');
% 
% k = pd.ParameterValues(1);  % 形状参数
% theta = pd.ParameterValues(2);  % 尺度参数
% 
% x = linspace(min(init_eta1), max(init_eta1), 100);
% y = pdf(pd, x);
% 
% figure;
% histogram(init_eta1, 'Normalization', 'pdf');  % 绘制直方图
% hold on;
% plot(x, y, 'r-', 'LineWidth', 2);  % 绘制拟合的 Gamma 分布
% hold off;
% title('Gamma Distribution Fit');
% xlabel('Data');
% ylabel('Probability Density');
