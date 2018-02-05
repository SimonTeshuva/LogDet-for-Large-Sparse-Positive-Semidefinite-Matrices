close all;
clear all;
clc;

a = 1e-3;
b = 1;
N = 15;
[t, w] = lgwt(N, a,b);

x = a:1e-3:b;
logx = log(x);

for i = 1:length(x)
logx_gauss(i) = (x(i)-1) * sum(w./(t*x(i) + (1-t)));
end

logx_gauss2 = (x-1) .* sum(w./(t.*x + (1-t)));


plot(x, logx, 'kd', x, logx_gauss, 'r-')
legend('log(x)', 'GQ log(x)')