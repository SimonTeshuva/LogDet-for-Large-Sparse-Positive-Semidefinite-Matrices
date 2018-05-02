close all;
clear all;
clc;

diag_dom_const = 1e-3;
samp = 20;
order = 14;
kappa = 5;
max_chol_size = 6e4;
max_chol_index = 5;

A_size_vect = [1e3, 5e3, 1e4, 2e4, 4e4, 1e5, 2e5, 4e5, 6e5, 8e5, 1e6];
A_sub_length = length(A_size_vect);

A_size_vect = A_size_vect(1:A_sub_length);
maxk = 11;

kneser_sizes = zeros(size(1,maxk));
for k=1:maxk
    kneser_sizes(k) = nchoosek(2*k+1,k);
end


logdet_exact_kneser = zeros(size(kneser_sizes));

rat_logdets_kneser = zeros(size(kneser_sizes));
rat_logdets_kneser_time = zeros(size(kneser_sizes));
cheb_logdets_kneser_time = zeros(size(kneser_sizes));
cheb_logdets_kneser = zeros(size(kneser_sizes));
chol_logdet_kneser = zeros(size(kneser_sizes));
chol_logdet_kneser_time = zeros(size(kneser_sizes));

rat_logdet_rand = zeros(size(A_size_vect));
rat_logdets_rand_time = zeros(size(A_size_vect));
cheb_logdets_rand = zeros(size(A_size_vect));
cheb_logdets_rand_time = zeros(size(A_size_vect));
chol_logdets_random = zeros(size(A_size_vect));
chol_logdets_random_time = zeros(size(A_size_vect));


block_size = min(20,samp);
counter=0;
while counter < length(A_size_vect)
    counter = counter+1;
    fileName = strcat('kneser (',num2str(2*counter+1),',',num2str(counter),').txt');
    [n,k,A_kneser] = generate_adjacency_matrix(fileName);
    fprintf('testing kneser graph matrix of for (n,k) = (%.f,%.f)\n', n,k);
    [evals,mults]=kneser_eigenvalues(n,k);

    
    lmin = -nchoosek(n-k-1,k);
    lmax = nchoosek(n-k,k);
    epsilon = (lmax-lmin)/(kappa+1);
    A_kneser_shift = diag(A_kneser)+speye(size(A_kneser))*(epsilon-lmin);
    
    fprintf('computing exact logdet of the kneser graph\n');
    logdet_exact_kneser(counter) = sum(mults.*log(evals));
    
    parfor i=1:2
        fprintf('starting parallel pool\n');
    end
    
    fprintf('computing logdet of the kneser graph using the rational function method\n');
    tic;
    rat_logdets_kneser(counter) = Rational(A_kneser_shift,order,samp,block_size);
    rat_logdets_kneser_time(counter) = toc;
    
    fprintf('computing logdet of the kneser graph using the chebyshev function method\n');
    tic;
    cheb_logdets_kneser(counter) = ChebLogDet(A_kneser_shift,samp,order,min(evals),max(evals));
    cheb_logdets_kneser_time(counter) = toc;
    if nchoosek(n,k)<=max_chol_size
        fprintf('computing logdet of the kneser graph using the cholesky decomposition method\n');
        tic
        amd_kneser = symamd(A_kneser_shift);
        chol_logdet_kneser(counter) = full(2*sum(log(diag(chol(A_kneser_shift(amd_kneser,amd_kneser))))));
        chol_logdet_kneser_time(counter) = toc;
    end
    
    clear A_kneser
    clear A_kneser_shift
    
    A_size = A_size_vect(counter);
    fprintf('testing random matrix of size %.f\n', A_size);
    A_rand = Sparse_Dataset(A_size, 2*floor(log10(A_size)),diag_dom_const);
    
    parfor i=1:2
        fprintf('starting parallel pool\n');
    end
    
    fprintf('computing logdet of the random matrix using the rational function method\n');
    tic;
    rat_logdet_rand(counter) = Rational(A_rand,order,samp,block_size);
    rat_logdets_rand_time(counter) = toc;
    
    fprintf('computing logdet of the random matrix using the chebyshev function method\n');
    tic;
    cheb_logdets_rand(counter) = ChebLogDet(A_rand,samp,order,max(max(A_rand)),diag_dom_const);
    cheb_logdets_rand_time(counter) = toc;
    
    if A_size<=max_chol_size
        fprintf('computing logdet of the random matrix using the choleskty decomposition method\n');
        tic
        amd_rand = symamd(A_rand);
        chol_logdets_random(counter) = full(2*sum(log(diag(chol(A_rand(amd_rand,amd_rand))))));
        chol_logdets_random_time(counter) = toc;
    end
    
    clear A_rand
    
end

% cheb_logdets_kneser is returning imaginary numbers


% 1a random (runtime)
figure(1)
plot(cheb_logdets_rand_time,log10(A_size_vect),rat_logdets_rand_time,log10(A_size_vect));
ylabel('log10(size)');
xlabel('time');
title('time vs matrix size for randomise matrix');
legend('Chebyshev','Rational');

% 1a kneser (runtime)
figure(2)
plot(cheb_logdets_kneser_time,log10(kneser_sizes),rat_logdets_kneser_time,log10(kneser_sizes));
ylabel('log10(size)');
xlabel('time');
title('time vs matrix size for kneser graph based matrix');
legend('Chebyshev','Rational');

% 1b random (relative accuracy vs chol)
rel_error_chol_cheb_rand = (chol_logdets_random - cheb_logdets_rand(1:length(chol_logdets_random)))./chol_logdets_random;
rel_error_chol_cheb_kneser = (chol_logdet_kneser - cheb_logdets_kneser(1:length(chol_logdet_kneser)))./chol_logdet_kneser;
rel_error_chol_rat_rand = (chol_logdets_random - rat_logdet_rand(1:length(chol_logdets_random)))./chol_logdets_random;
rel_error_chol_rat_kneser = (chol_logdet_kneser - rat_logdets_kneser(1:length(chol_logdet_kneser)))./chol_logdet_kneser;
A_size_sub = A_size_vect(1:length(rel_error_chol_cheb_rand));

figure (3)
plot(log10(A_size_sub), rel_error_chol_cheb_rand, log10(A_size_sub),rel_error_chol_rat_rand);
xlabel('log10(size)');
ylabel('relative accuracy');
title('relative accuracy (compared to cholesky decomposition) vs matrix size for randomised matrix');
legend('Chebyshev','Rational');

% 1b kneser (relative accuracy vs chol)
figure (4)
plot(log10(A_size_sub), rel_error_chol_cheb_kneser(1:length(A_size_sub)), log10(A_size_sub),rel_error_chol_rat_kneser(1:length(A_size_sub)));
xlabel('log10(size)');
ylabel('relative accuracy');
title('relative accuracy (compared to cholesky decomposition) vs matrix size for matrix based on kneser graph');
legend('Chebyshev','Rational');

% 1c random (runtime vs chol)
figure(5)
plot(cheb_logdets_rand_time(1:max_chol_index),log10(1:max_chol_index),rat_logdets_rand_time(1:max_chol_index),log10(1:max_chol_index),chol_logdets_random_time(1:max_chol_index), log10(1:max_chol_index));
ylabel('log10(size)');
xlabel('time');
title('time vs matrix size for random matrix');
legend('Chebyshev','Rational','Cholesky');

% 1c kneser (runtime vs chol)
figure(6)
plot(cheb_logdets_kneser_time(1:length(A_size_sub)),log10(A_size_sub),rat_logdets_kneser_time(1:length(A_size_sub)),log10(A_size_sub),chol_logdet_kneser_time(1:length(A_size_sub)),log10(A_size_sub));
ylabel('log10(size)');
xlabel('time');
title('time vs matrix size for kneser graph based matrix');
legend('Chebyshev','Rational','Cholesky');


% 1d kneser (relative accuracy vs kneser)
rel_error_cheb_kneser = (logdet_exact_kneser - cheb_logdets_kneser)./logdet_exact_kneser;
rel_error_rat_kneser = (logdet_exact_kneser - rat_logdets_kneser)./logdet_exact_kneser;

figure (7)
plot(log10(A_size_vect), rel_error_cheb_kneser(1:length(A_size_vect)), log10(A_size_vect),rel_error_rat_kneser(1:length(A_size_vect)))
xlabel('log10(size)');
ylabel('relative accuracy');
title('relative accuracy vs matrix size for matrix based on kneser graph');
legend('Chebyshev','Rational');

save('results.mat');

function [eigenvalues, multiplicities] = kneser_eigenvalues(n,k)
evals = [];
mults = [];
for i = 0:k
    eval = ((-1)^i)*nchoosek(n-k-i, k-i);
    if i == 0
        mult = 1;
    else
        mult = nchoosek(n,i) - nchoosek(n, i-1);
    end
    evals = [evals eval];
    mults = [mults mult];
end
eigenvalues = evals;
multiplicities = mults;
end