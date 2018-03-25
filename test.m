function test_data = test()
close all;
clear all;
clc;

disp('starting');

dataset_size = 1e3;
epr = 10;
diag_dom_const = 1e-3;
tic
dataset = Sparse_Dataset(dataset_size, epr, diag_dom_const);
dataset_time = toc;
disp('dataset done')
n = 14;
Nm = 100;

M = 20;
a = 1e-3;
b = 1;
k = 5;

logdet_cheb = [];
logdet_rat = [];
logdet_rat_re = [];
logdet_exact = [];
for epr = 5:5:20
    for y = 1:3:7
        dataset_size = y*10^3;
        dataset = Sparse_Dataset(dataset_size, epr, diag_dom_const);
        logdet_exact = [logdet_exact full(2*sum(log(diag(chol(dataset)))))];
        logdet_rat = [logdet_rat Parallel_Rational(dataset, M, Nm)];
        logdet_rat_re = [logdet_rat_re Parallel_Rational_Reordered(dataset, M, Nm, k)];
        logdet_cheb = [logdet_cheb Approx_Algorithm(dataset, Nm, n, diag_dom_const)];
    end
end

test_data = [logdet_exact; logdet_rat; logdet_rat_re; logdet_cheb];
end