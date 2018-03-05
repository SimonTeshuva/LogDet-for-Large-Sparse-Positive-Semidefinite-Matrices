close all;
clear all;
clc;
% https://au.mathworks.com/help/matlab/ref/sprandsym.html
%need a far better method. this bricks my computer at somewhere between 1E4
%and 2E4

datset_size = 1e3;
epr = 10;
diag_dom_const = 1e-3;
dataset = Sparse_Dataset(datset_size, epr, diag_dom_const);
n = 5;
m = 10;
%  logdet = Approx_Algorithm(dataset, m, n, diag_dom_const);
% fprintf('logdet(A) = %.5d + %.5di\n', real(logdet), imag(logdet));

R = chol(dataset);
exact_val = 2*sum(log(diag(R)));
cheb_approx = Approx_Algorithm(dataset, m, n, diag_dom_const);

M = 50;
a = 1e-3;
b = 1;
N = 50;
rational_approx = New_Algorithm(dataset, M, N);

% for N=5:5:50
%     for M = 5:5:50
%         logdet_new = New_Algorithm(dataset, M, N);
%         fprintf('M = %f, N=%f, logdet = %f\n', M,N,logdet_new);
%     end
% end

% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory