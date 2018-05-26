%% Rational Function Approximator for Logdet
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% This function takes in a large, sparse, positive definite matrix, 
% an order for the Rational Function Approximator (M), the number of 
% samples to be averaged (N), and the size of the block to be computed in each 
% iteration (k) and estimates the logdet of the matrix. 

function logdet = Rational(dataset, M, N, k)
% optimal values for N and k were found to be 20 and 20 on my computer

B = dataset;
dataset_size = size(B);
n = dataset_size(1);
I = speye(n,n);
[t, w] = lgwt(M, 0, 1); % generate nodes and weights for the rational function
% approximator
V = ((rand(N,n)<.5)*2 - 1)'; % N randmacher vectors;

% derivation for the formula to compute the logdet
% logdet=(1/N)*SUM(1:N){v(i)'*(B-I)*SUM(1:M){wj*(tj*B+(1-tj)*I)^(-1)*v(i)}}
% logdet
% =SUM(1:M){(wj/N)*SUM(1:N){v(i)'*(B-I)*pcg(tj*B+(1-tj)*I),v(i)}}
% =SUM(1:M){(wj/N)*SUM(1:N){v(i:(i+k))'*(B-I)*pcg(tj*B+(1-tj)*I),v(i:(i+k))}}

% change N a little bit so that k will be a factor of N
fake_N = k*floor(N/k); 


logdet_sum = 0; % this will store the sum of the samples for the logdet
parfor j = 1:M % run in parallel because each iteration is independant of the
    % others. this is because the rational function can be broken up into a
    % set of partial fractions
    tj = t(j);
    wj = w(j);
    total = 0;
    for i = 1:k:fake_N % solve fake_N many linear systems. solve k at a time
        % logdet=SUM(1:M){(wj/N)*SUM(1:N){v(i:(i+k))'*(B-I)*pcg(tj*B+(1-tj)*I),v(i:(i+k))}}

        v = V(:, i:(i+k-1)); % take k of the rademacher vectors

        % break down the formula for estimating the logdet into a number of
        % steps
        A = tj*B + (1-tj)*I; 
        afun = @(v) reshape(A*(reshape(v, n, k)),n*k, 1); % reshaping A and 
        % v so that they can be used as inputs in pcg. using a function
        % handle to do so.
        pcg_ret = reshape(pcg_quiet(afun, v(:)), n, k);
        BIpcg = (B-I)*pcg_ret;
        
        LS_solved = wj*v(:)'*BIpcg(:);
        total = total + LS_solved;
    end
    logdet_sum = logdet_sum + total;
    
    fprintf('%.f\n', 100*j/M);
end
logdet = logdet_sum/fake_N; % take the average of the fake_N many samples
end