% 计算lambda(t)系数初始值
%%
function  [Lambda_hat, lambda_hat, delta, gamma, ETA, a, beta] = initial_Lambda(Y, n, model, Ct, init_delta, init_gamma, init_beta, Cox_z)

for iter = 1:1e3

for i = 1:n
    Y{i} = [Y{i}(1), diff(Y{i})];
end

for i = 1:n
    for j = 1:length(Ct)
        y(i,j) = Y{i}(j);
    end
end
y_mean = mean(y);
for j = 1:length(Ct)
    a = sum(y);
    init_lambda_hat(j) = y(j);
    % init_lambda_hat(j) = 0.2;
end

fprintf('初始迭代 \n')

for i = 1:n
    total = 0;
    for j = 1:length(Ct)
        total = total + (Y{i}(j) - init_lambda_hat(j)*exp(init_beta*Cox_z(j+1+(17*(i-1))))).^2./(2*Y{i}(j));
        % disp(total)
    end
    log_eta(i) = psi(init_delta + length(Ct)/2) - log(init_gamma + total);
    eta(i) = exp(log_eta(i));
    disp(eta(i))
end
% 
% eta = 7 + 16*rand(1, 15);

for j = 1:length(Ct)
    options = optimoptions('fsolve', 'Display', 'none'); % 显示迭代过程
    initial_guess = [init_lambda_hat, init_beta];
    solution = fsolve(@(x) nonlinear_equations1(x, eta, Y, n, Ct, Cox_z), initial_guess, options);
    for j = 1:length(Ct)
        lambda_hat(j) = solution(j);
    end
    beta = solution(length(Ct)+1);
end

lb = [0.01, 0.1];  % gamma 和 delta 不能为负，防止过小
ub = [10e5, 10e6];    % 防止过大导致 eta_i 退化

% 设定 fmincon 选项
options = optimoptions('fsolve', 'Display', 'none');
x0 = [init_gamma, init_delta];
% 用 fmincon 进行优化求解
% [x, ~] = fmincon(@(x) myNonlinearSystem(x, eta), x0, [], [], [], [], lb, ub, [], options);
x = fsolve(@(x) myNonlinearSystem(x, eta), x0, options);

% 输出结果
gamma = x(1);
delta = x(2);

if gamma > 1e6
    init_delta = delta;
    continue
else
    break
end


end

for j = 1:length(Y{1})
    if j == 1
        Lambda_hat(j) = lambda_hat(j);
    else
        Lambda_hat(j) = lambda_hat(j) + Lambda_hat(j-1);
    end
end

ETA = eta;

Ct_2 = [0,Ct];
switch model
    case 'linear'
        a = sum(Ct_2 .* [0,Lambda_hat]) / sum(Ct_2 .^2);
    case 'quadratic'
        a = sum(Ct_2 .^2./10 .* [0,Lambda_hat]) / sum(Ct_2.^4/100);
    case 'sqrt'
        a = sum(sqrt(10 .* Ct_2 ).* [0,Lambda_hat]) / sum(10 .* Ct_2);
    otherwise
        error('Unknown model type');
end

end

%%
function y = nonlinear_equations1(x, eta, Y, n, Ct, Cox_z)

for j = 1:length(Ct)
    sum_term(j) = 0;
    for i = 1:n
        sum_term(j) = sum_term(j) + 1 + x(j)*eta(i)*exp(x(length(Ct)+1)*Cox_z(j+1+(17*(i-1))))*(1 - x(j)./Y{i}(j)*exp(x(length(Ct)+1)*Cox_z(j+1+(17*(i-1)))));
    end
end

sum_term1 = 0;
for i = 1:n
    sum_term2 = 0;
    for j = 1:length(Ct)
        sum_term2 = sum_term2 + Cox_z(j+1+(17*(i-1)),1) + Cox_z(j+1+(17*(i-1)),1)*x(j)*eta(i)*exp(x(length(Ct)+1)*Cox_z(j+1+(17*(i-1))))*(1-x(j)./Y{i}(j)*exp(x(length(Ct)+1)*Cox_z(j+1+(17*(i-1)))));
    end
    sum_term1 = sum_term1 + sum_term2;
end


for j = 1:length(Ct)
    y(j) = sum_term(j);
end
y(length(Ct)+1) = sum_term1;

end
%%

function F = myNonlinearSystem(x, eta)
    gamma = x(1); % 第一个变量 γ
    delta = x(2); % 第二个变量 δ
    n = length(eta);

    % 第一方程
    % F(1) = sum(log(gamma) + log(eta) - (log(delta-1./(2*delta)-1./(12*delta.^2))));
    F(1) = sum(log(gamma) + log(eta) - (psi(delta)));
    % 第二方程
    F(2) = sum(delta / gamma - eta);

    % obj = F1 + F2 + 0.01*(delta/gamma).^2;
end

