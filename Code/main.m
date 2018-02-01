close all;
clear all;
clc;
% https://au.mathworks.com/help/matlab/ref/sprandsym.html
%need a far better method. this bricks my computer at somewhere between 1E4
%and 2E4

dataset = Sparse_Dataset(1e3, 10);
logdet = Approx_Algorithm(dataset, 10, 10);
fprintf('logdet(A) = %.5d + %.5di\n', real(logdet), imag(logdet));

exact_val = log(sum(eig(dataset)))

% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory