f = csvread("noisy_signal.csv");
omega_val = [1e-4, 5e-4, 1e-3];
row_num = 999;
col_num = 1000;
f_step = 1/col_num;
f_range = f_step:f_step:1;
L = zeros(row_num, col_num);
for i = 1:row_num
    L(i,i) = row_num;
    L(i,i+1) = -row_num;
end
for i = 1:length(omega_val)
    x = (eye(col_num) + omega_val(i)*(L'*L))\f;
    figure;
    plot(f_range,f,'k.',f_range,x,'r.');
    xlabel('$x$','Interpreter','latex'),ylabel('$s(x)$','Interpreter','latex');
    title(['Denoised signal with \omega=',num2str(omega_val(i))]);
    legend('Noisy signal','Denoised Signal');
end




