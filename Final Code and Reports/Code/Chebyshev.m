%% Chebyshev Function Approximation for the Logdet
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% This function uses a chebyshev function approximator to estimate the
% logdet of a large, sparse, positive definite matrix. it takes as input; a
% matrix to be evaluated (dataset), a number of samples to be averaged (m),
% an order for the chebyshev function (n).

% it uses gershgorin circle theorem to get an estimate for the maximum and
% minimum eigenvalues of the matrix, which are then used to normalize the
% matrix so that its eigenvalues lie in the range [-1,1]. A chebyshev
% function approximator is then generated using the CHEBFUN library.

% at this point, the algorithm devised by I. Han, D. Malioutov, J. Shin in
% “Large-scale Log-determinant Computation through Stochastic Chebyshev Expansions”, 
% is implemented to estimate the logdet, and the result is then
% de-nomralised and returned as output. 
function logdet = Chebyshev(dataset, m, n)
format long
dataset_size = size(dataset);


B = dataset;
I = speye(size(B));
A = I - B;

% normalizin the matrix A to have eigenvalues in the range that Chebyshev
% Functions are accurate [-1,1]. 

% sigma_max = smax, sigma_min = smin
% [smin .... smax]
%  v  v  v v v 
% [delta ... 1-delta]
% --> alpha*smin = delta, alpha*smax = 1-delta
% solve:
% alpha = 1/(smax+smin)
% delta = alpha*smin

smax = max(max(dataset)); % using Gershgorin Circle Theorem
smin = min(2*diag(dataset)'-sum(abs(A))); % using Gershgorin Circle Theorem

alpha = 1/(smax+smin);
delta = alpha*smin;

A_new = alpha*A;

% using http://www.chebfun.org/ to generate chebyshev polynomial
% coeficcients
x = chebfun('x');
p = log(1- ((1-2*delta)*x+1)/2); % using this d we get a complex chebyshev fn
c = chebpoly(p,n); % this is not the funcion you are looking for
% plot c and p to see if we make sense
% investigate chebfun
% normalisation of chebyshev coefficients

G = 0; %will store the normalized estimate for the logdet

for i = 1:m % compute m sample estimates for the logdet
    v = ((rand(1,dataset_size(1))<.5)*2 - 1)'; % randmacher vector;
    u = c(1)*v;
    if n>1 % for each sample, estimate the logdet using a chebyshev 
        % function approximator of order n
        w0 = v;
        w1 = A_new*v;
        u = u + c(2)*w1;
        for j = 3:(n)
            w2 = 2*A_new*w1 - w0; % derived from the chebyshev function
            u = u + c(j)*w2;
            w0 = w1;
            w1 = w2;
        end
    end
    fprintf('%.f\n', 100*i/m);
    G = G + v'*u/m;
end

% de-normalizing the result

% logdet(A_new)=logdet(alpha*A_old)
%              =log(det(alpha*I)*det(A_old))
%              =size*log(alpha)+logdet(A_old)
% logdet(A_old)=logdet(A_new)-size*log(alpha)
logdet = G-dataset_size(1)*log(alpha);
end