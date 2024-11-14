f = readmatrix("noisy_signal.csv");
row_num = 999;
col_num = 1000;
f_step = 1/col_num;
f_range = f_step:f_step:1;
L = zeros(row_num, col_num);
for i = 1:row_num
    L(i,i) = row_num;
    L(i,i+1) = -row_num;
end
lambda = 5e-4;
x0 = f;
y0 = L*f;
z = ones(row_num,1);
u_old = zeros(row_num+col_num,1);
u_new = [x0;y0];
function [u_norm] = u_approx(u_new)
row_num = 999;
col_num = 1000;
omega = 75;
gamma = 1000;
u_norm = 0;
for i = (1 + col_num):(row_num + col_num)
    if (abs(u_new(i,1)) <= 1/gamma)
        u_norm = u_norm + omega*0.5*gamma*u_new(i,1)^2;
    else
        u_norm = u_norm + omega*(abs(u_new(i,1))-1/(2*gamma));
    end
end
end

function [u_grad] = grad_u_approx(u_new)
omega = 75;
gamma = 1e3;
row_num = 999;
col_num = 1000;
u_grad = zeros(row_num + col_num,1);
for i = (col_num + 1):(col_num + row_num)
    if (abs(u_new(i,1)) <= 1/gamma)
        u_grad(i,1) = omega*gamma*u_new(i,1);
    elseif (u_new(i,1) < -1/gamma)
        u_grad(i,1) = -omega;
    else
        u_grad(i,1) = omega;
    end
end
end

A_bar = [eye(col_num),zeros(col_num,row_num)];
L_bar = [-L,eye(row_num)];
func_f = @(u_new) 0.5*norm(A_bar*u_new - f)^2 + 0.5*lambda*norm(L_bar*u_new-z)^2+u_approx(u_new);
func_f_grad = @(u_new) A_bar'*(A_bar*u_new-f) + lambda*L_bar'*(L_bar*u_new - z) + grad_u_approx(u_new);


alpha = 0.1;
beta = 0.5;
iter_num = 0;
step = 1;
tic;
while iter_num < 500 && norm(u_new-u_old)/norm(u_new) >= 1e-4
    iter_num = iter_num + 1;
    while func_f(u_new)-func_f(u_new-step*func_f_grad(u_new)) < alpha*step*norm(func_f_grad(u_new))^2
        step = step*beta;
    end
    u_old = u_new;
    u_new = u_new - step*func_f_grad(u_new);
end

z = z + L*u_old(1:col_num)-u_old(col_num+1:col_num+row_num);
while norm(u_new(1:col_num)-u_old(1:col_num))/norm(u_new(1:col_num)) >= 1e-5
    iter_num = iter_num + 1;
    while func_f(u_new)-func_f(u_new-step*func_f_grad(u_new)) < alpha*step*norm(func_f_grad(u_new))^2
        step = step*beta;
    end
    u_old = u_new;
    u_new = u_new - step*func_f_grad(u_new);
end
time = toc;
figure;
plot(f_range, f, 'k.',f_range, u_old(1:col_num,1), 'r.');
xlabel('$x$','Interpreter','latex'),ylabel('$s(x)$','Interpreter','latex');
title('Denoised Signal vs Noisy signal');
legend('Noisy signal','Denoised Signal');
fprintf('The number of iterations is %d, execution time is %.4f seconds\n',iter_num,time);

