function [Lambda_hat, lambda_hat, init_a, eta, delta, gamma, beta, Lambda_new] = EM_algorithm_hat1(Y, TNi, init_a, tol, max_iter, model, n, eta1, Ct, init_lambda_hat, init_lambda, init_delta, init_gamma, init_beta, Cox_z)

init_a_tol = zeros(max_iter+1, 1);
init_a_tol(1) = init_a;
Eta{1} = eta1;
Delta(1) = init_delta;
Gamma(1) = init_gamma;
Beta(1) = init_beta;
Ct_2 = [0, Ct];
Lambda_new = cell(max_iter+1, 1);
lambda_new = cell(max_iter+1, 1);
lambda_new{1} = init_lambda_hat;
Lambda_new{1} = init_lambda;

for iter = 1:max_iter

    % E-step: 计算期望

    [E_yij, E_inv_yij] = E_step(Y, init_a_tol(1), TNi, model, Eta(iter), n);
    for i = 1:n
        for j = 1:length(E_yij{i})
            E_Yij(i,j) = E_yij{i}(j);
        end
    end

    % M-step: 最大化期望
    fprintf('EM 迭代次数为%d \n', iter)

    [Lambda_new{iter+1}, lambda_new{iter+1}, Eta{iter+1}, Delta(iter+1), Gamma(iter+1), Beta(iter+1)] = M_step(E_yij, E_inv_yij, n, Eta{iter}, Delta(iter), Gamma(iter), Y, Ct, Beta(iter), lambda_new{iter}, Cox_z);
    
    

    switch model
        case 'linear'
            init_a_tol(iter+1) = sum(Ct_2 .* [0,Lambda_new{iter+1}]) / sum(Ct_2 .^2);
        case 'quadratic'
            init_a_tol(iter+1) = sum(Ct_2 .^2./10 .* [0,Lambda_new{iter+1}]) / sum(Ct_2.^4/100);
        case 'sqrt'
            init_a_tol(iter+1) = sum(sqrt(10 .* Ct_2 ).* [0,Lambda_new{iter+1}]) / sum(10 .* Ct_2);
        otherwise
            error('Unknown model type');
    end

    if abs(sum((Lambda_new{iter+1}-Lambda_new{iter}).^2 * 0.92)) < tol
        init_a = init_a_tol(iter+1);
        eta = Eta{iter+1};
        delta = Delta(iter+1);
        gamma = Gamma(iter+1);
        beta = Beta(iter+1);
        Lambda_hat = Lambda_new{iter+1};
        lambda_hat = lambda_new{iter+1};
        break
    else
        continue
    end

end

end

