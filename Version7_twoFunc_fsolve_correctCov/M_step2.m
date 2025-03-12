function [Lambda_hat, lambda_hat, ETA, delta, gamma] = M_step2(E_yij, E_inv_yij, n, init_eta, init_delta, init_gamma, Y, Ct)
for i = 1:n
    Y{i} = [Y{i}(1), diff(Y{i})];
end

x0 = [1e-8, 1e2];
for j = 1:length(Ct)
    lambda_hat(j) = fzero(@(lambda_j) equation18_lambdaj(init_eta, n, Y, j, lambda_j), x0);
end

for i = 1:n
    total = 0;
    for j = 1:length(Ct)
        total = sum((Y{i} - lambda_hat).^2./(2*Y{i}));
    end
    log_eta(i) = psi(init_delta + length(Ct)/2) - log(init_gamma + total);
    eta(i) = exp(log_eta(i));
end

% z1 = sym('z1');
% z2 = sym('z2');
% eq1 = n*log(z2)+sum(log(eta))-n*psi(z1)==0;
% eq2 = n*z1./z2-sum(eta)==0;
% [delta,gamma] = solve([eq1,eq2],[z1,z2]);


number = 100000;
Delta = zeros(number,1);
Gamma = zeros(number,1);
Delta(1) = init_delta;
Gamma(1) = init_gamma;
for j = 2:number
    Delta(j) = fzero(@(delta)fun_delta(delta, Gamma(j-1),eta,n), Delta(j-1));
    Gamma(j) = fzero(@(gamma)fun_gamma(Delta(j),gamma,eta,n), Gamma(j-1));
    if abs(Delta(j) - Delta(j-1)) < 1
        delta = Delta(j);
        gamma = Gamma(j);
    else
        continue
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

end
%%
function y = equation18_lambdaj(eta, n, Y, j, lambda_j)
sum_term = 0;
for i = 1:n
    sum_term = sum_term + eta(i)-eta(i) * lambda_j ./Y{i}(j);
end
y = n./ lambda_j + sum_term;
end


function y = fun_delta(delta, gamma, eta, n)
y = n*log(gamma) + sum(log(eta)) - n*psi(delta);
end

function y = fun_gamma(delta, gamma, eta, n)
y= n*delta./gamma - sum(eta);

end
