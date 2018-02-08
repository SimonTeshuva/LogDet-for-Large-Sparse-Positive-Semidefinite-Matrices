close all;
clear all;
clc;

tmin = 0;
tmax = 1;
N = 25; % play around with
%  MATLAB log? what does it do
[t, w] = lgwt(N,tmin,tmax);

a = 1e-3;
b = 1/a;
x = a:1e-3:b;
logx = log(x);

for i = 1:length(x)
logx_gauss(i) = (x(i)-1) * sum(w./(t*x(i) + (1-t)));
end

logx_gauss2 = (x-1) .* sum(w./(t.*x + (1-t)));


% plot(x, logx, 'kd', x, logx_gauss, 'r-')
error = logx_gauss - logx;
loglog(x, abs(error));
legend('log(x)', 'GQ log(x)')