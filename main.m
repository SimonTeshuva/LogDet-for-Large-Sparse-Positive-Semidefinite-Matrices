close all;
clear all;
clc;
% https://au.mathworks.com/help/matlab/ref/sprandsym.html
%need a far better method. this bricks my computer at somewhere between 1E4
%and 2E4

disp('starting');

dataset_size = 1e4;
epr = 10;
diag_dom_const = 1e-3;
tic
dataset = Sparse_Dataset(dataset_size, epr, diag_dom_const);
dataset_time = toc;
disp('dataset done')
n = 14;
m = 15;
%  logdet = Approx_Algorithm(dataset, m, n, diag_dom_const);
% fprintf('logdet(A) = %.5d + %.5di\n', real(logdet), imag(logdet));

if dataset_size < 6e3
    tic
    R = chol(dataset);
    exact_val = full(2*sum(log(diag(R))));
    exact_time = toc;
    disp('exact done');
end
tic
cheb_approx = Approx_Algorithm(dataset, m, n, diag_dom_const);
cheb_time = toc;
disp('cheb done');

M = 20;
a = 1e-3;
b = 1;
N = 100;
% tic
% rational_approx = New_Algorithm(dataset, M, N);
% rational_time = toc;
% disp('rat done');

tic
prational_approx = Parallel_Rational(dataset, M, N);
prational_time = toc;
disp('prat done');


clc
if dataset_size < 6e4
fprintf('result from %s method is %.3f, and took %.3f milliseconds \n', 'exact', exact_val, 1000*exact_time);
end
fprintf('result from %s method is %.3f, and took %.3f milliseconds \n', 'chebyshev', cheb_approx, 1000*cheb_time);
% fprintf('result from %s method is %.3f, and took %.3f milliseconds \n', 'rational', rational_approx, 1000*rational_time);
fprintf('result from %s method is %.3f, and took %.3f milliseconds \n', 'parallel rational', prational_approx, 1000*prational_time);

% for N=5:5:50
%     for M = 5:5:50
%         logdet_new = New_Algorithm(dataset, M, N);
%         fprintf('M = %f, N=%f, logdet = %f\n', M,N,logdet_new);
%     end
% end

% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory