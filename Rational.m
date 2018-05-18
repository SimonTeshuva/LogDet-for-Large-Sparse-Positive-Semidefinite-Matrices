%% Rational Function Approximator for Logdet
% Author: Simon Tesuva
% Date Last Modified: 5/4/2018

% This function takes in a large, sparse, positive definite matrix, 
% an order for the Rational Function Approximator, the number of 
% samples desired, and the size of the block to be computed in each 
% iteration and computes the logdet of the matrix. 

function logdet = Rational(dataset, M, N, k)
% optimal N is 20 done in 1 block

B = dataset;
dataset_size = size(B);
n = dataset_size(1);
I = speye(n,n);
[t, w] = lgwt(M, 0, 1); % generate nodes and weights
V = ((rand(N,n)<.5)*2 - 1)'; % N randmacher vectors;

% derivation for the formula to compute the logdet
% logdet=(1/N)*SUM(1:N){v(i)'*(B-I)*SUM(1:M){wj*(tj*B+(1-tj)*I)^(-1)*v(i)}}
% logdet
% =SUM(1:M){(wj/N)*SUM(1:N){v(i)'*(B-I)*pcg(tj*B+(1-tj)*I),v(i)}}
% =SUM(1:M){(wj/N)*SUM(1:N){v(i:(i+k))'*(B-I)*pcg(tj*B+(1-tj)*I),v(i:(i+k))}}
fake_N = k*floor(N/k); % if N is not divisable by k, lower N such that N is
% divisable by k

gtotal = 0;
parfor j = 1:M
    tj = t(j);
    wj = w(j);
    total = 0;
    for i = 1:k:fake_N
        v = V(:, i:(i+k-1));
        A = tj*B + (1-tj)*I;
        afun = @(v) reshape(A*(reshape(v, n, k)),n*k, 1); % reshaping A and 
        % v so that they can be used as inputs in pcg. using a function
        % handle to do so.
        pcg_ret = reshape(pcg_quiet(afun, v(:)), n, k);
        BIpcg = (B-I)*pcg_ret;
        
        LS_solved = wj*v(:)'*BIpcg(:);
        total = total + LS_solved;
    end
    gtotal = gtotal + total;
    
    fprintf('%.f\n', 100*j/M);
end
logdet = gtotal/fake_N;
end