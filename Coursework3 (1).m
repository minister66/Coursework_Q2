f = csvread("noisy_signal.csv");
iter = 0;
lambda = 1e-2;
omega = 2e-3;
row_num = 999;
col_num = 1000;
L = zeros(row_num, col_num);
for i = 1:row_num
    L(i,i) = 999;
    L(i,i+1) = -999;
end
x_val_old = f;
y_val = L*x_val_old;
z_val = ones(row_num,1);
tic;
x_val_new = (eye(col_num) + lambda*(L'*L))\(f + lambda*L'*(y_val-z_val));
while norm(x_val_new - x_val_old)/norm(x_val_new) >= 1e-5
    for i = 1:row_num
        t = L(i,i)*x_val_new(i,1) + L(i,i+1)*x_val_new(i+1,1) + z_val(i,1);
        eta = omega/lambda;
        y_val(i,1) = (t/abs(t))*max(abs(t)-eta, 0);
    end
    z_val = z_val + L*x_val_new - y_val;
    iter = iter + 1;
    x_val_old = x_val_new;
    x_val_new = (eye(col_num) + lambda*(L'*L))\(f + lambda*L'*(y_val-z_val));
end
times = toc;
fprintf('The number of iterations is %d, and the execution time is %.4f seconds', iter, times);

