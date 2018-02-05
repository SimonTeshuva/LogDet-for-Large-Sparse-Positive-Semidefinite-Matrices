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
n = 10;
m = 10; 
logdet = Approx_Algorithm(dataset, m, n);

fprintf('logdet(A) = %.5d + %.5di\n', real(logdet), imag(logdet));

exact_val = log(sum(eig(dataset)))

% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory