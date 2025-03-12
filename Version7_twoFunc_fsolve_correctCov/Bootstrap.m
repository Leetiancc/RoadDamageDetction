function [A1, A1_Lambda_hat, A1_lambda_hat, ETA, Delta, Gamma] = Bootstrap(n, Ct, init_a, eta, model, tol, max_iter, lambda_new, B, TNi, Y, Lambda_new, init_delta, init_gamma)
yij_1_hat = zeros();
yij_2_hat = zeros();
E_yij_half_1 = cell(n,1);
E_yij_half_2 = cell(n,1);
P_1 = cell(n,1);
P_2 = cell(n,1);
delta = cell(n,1);
E_yij_Half = cell(n,1);
A1 = zeros();


for k = 1:B
    switch model
        case 'linear'
            lambda = @(t) init_a * t; % lambda 是 t 的函数
        case 'quadratic'
            lambda = @(t) init_a * (t.^2 / 10); % lambda 是二次函数
        case 'sqrt'
            lambda = @(t) init_a * (10*t)^0.5; % lambda 是平方根函数
        otherwise
            warning('Unknown model type');
    end

    % 生成样本
    for i = 1:n
        delta{i} = chi2rnd(1, 1, length(TNi{i}));
    end

    for i = 1:n
        for j = 1:length(TNi{i})
            if j == 1
                yij_1_hat(i,j) = 0.5 * (2*lambda(TNi{i}(j))+delta{i}(j)/eta(i) - sqrt(4*lambda(TNi{i}(j))*delta{i}(j)/(eta(i))+delta{i}(j)^2/(eta(i)^2)));
                P_1{i}(j) = lambda(TNi{i}(j)) / (lambda(TNi{i}(j)) + yij_1_hat(i,j));

                yij_2_hat(i,j) = 0.5 * (2*lambda(TNi{i}(j))+delta{i}(j)/eta(i) + sqrt(4*lambda(TNi{i}(j))*delta{i}(j)/(eta(i))+delta{i}(j)^2/(eta(i)^2)));
                P_2{i}(j)  = yij_1_hat(i,j) / (lambda(TNi{i}(j)) + yij_1_hat(i,j));
                
                r = rand;
                if r < P_1{i}(j)
                    E_yij_half_1{i}(j) = yij_1_hat(i,j);
                    E_yij_half_2{i}(j) = yij_2_hat(i,j);

                else
                    E_yij_half_1{i}(j) = yij_2_hat(i,j);
                    E_yij_half_2{i}(j) = yij_2_hat(i,j);

                end

            else
                yij_1_hat(i,j) = 0.5 * (2*(lambda(TNi{i}(j))-lambda(TNi{i}(j-1)))+delta{i}(j)/eta(i)-sqrt(4*(lambda(TNi{i}(j))-lambda(TNi{i}(j-1)))*delta{i}(j)/(eta(i))+delta{i}(j)^2/(eta(i)^2)));
                P_1{i}(j)  = (lambda(TNi{i}(j))-lambda(TNi{i}(j-1))) / ((lambda(TNi{i}(j))-lambda(TNi{i}(j-1))) + yij_1_hat(i,j));

                yij_2_hat(i,j) = 0.5 * (2*(lambda(TNi{i}(j))-lambda(TNi{i}(j-1)))+delta{i}(j)/eta(i)+sqrt(4*(lambda(TNi{i}(j))-lambda(TNi{i}(j-1)))*delta{i}(j)/(eta(i))+delta{i}(j)^2/(eta(i)^2)));
                P_2{i}(j)  = yij_1_hat(i,j) / ((lambda(TNi{i}(j))-lambda(TNi{i}(j-1))) + yij_1_hat(i,j));

                r = rand;
                if  r < P_1{i}(j)
                    E_yij_half_1{i}(j) = yij_1_hat(i,j);
                    E_yij_half_2{i}(j) = yij_2_hat(i,j);

                else
                    E_yij_half_1{i}(j) = yij_2_hat(i,j);
                    E_yij_half_2{i}(j) = yij_1_hat(i,j);

                end
            end
        end
    end  

    for i = 1:n
        switch model
            case 'linear'
                E_yij_Half{i} = E_yij_half_1{i};
            case 'quadratic'
                E_yij_Half{i} = E_yij_half_1{i};
            case 'sqrt'
                E_yij_Half{i} = E_yij_half_1{i};
            otherwise
                warning('Unknown model type');
        end
    end

    [Lambda_hat, lambda_hat, eta, delta1, gamma] = EM_algorithm_hat2(Y, TNi, init_a, tol, max_iter, model, n, eta, Ct, lambda_new, E_yij_Half, Lambda_new, init_delta, init_gamma);

    Ct_2 = [0, Ct];
    switch model
        case 'linear'
            a = sum(Ct_2 .* [0,Lambda_hat]) / sum(Ct_2 .^2);
        case 'quadratic'
            a = sum(Ct_2 .^2./10 .* [0,Lambda_hat]) / sum(Ct_2.^4/100);
        case 'sqrt'
            a = sum(sqrt(10 .* Ct_2 ).* [0,Lambda_hat]) / sum(10 .* Ct_2);
        otherwise
            warning('Unknown model type');
    end

    A1(k) = a;
    ETA{k} = eta;
    Delta(k) = delta1;
    Gamma(k) = gamma;
    A1_Lambda_hat{k} = Lambda_hat;
    A1_lambda_hat{k} = lambda_hat;

    fprintf('bootstrap迭代次数为: %d\n ', k)

end

end





