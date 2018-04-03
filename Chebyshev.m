function logdet = Chebyshev(dataset, m, n, diag_dom_const)
format long
dataset_size = size(dataset);


B = dataset;
I = speye(size(dataset));
A = I - B;

% sigma_max = smax, sigma_min = smin
% [smin .... smax]
%   v  v  v v v 
% [delta ... 1-delta]
% --> alpha*smin = delta, alpha*smax = 1-delta
% solve:
% alpha = 1/(smax+smin)
% delta = alpha*smin

% smax = max(max(diag(diag(A)))); 
% smin = min(min(diag(diag(A))));
smax = 2*max(max(dataset)); % max(sub(abs(A)))
smin = diag_dom_const; % min(2*diag(A) - sum(abs(A)))

alpha = 1/(smax+smin);
delta = alpha*smin;

A_new = alpha*A;

G = 0;
% using http://www.chebfun.org/ to generate chebyshev polynomial
% coeficcients

x = chebfun('x');
p = log(1- ((1-2*delta)*x+1)/2); % using this d we get a complex chebyshev fn
c = chebpoly(p,n); % this is not the funcion you are looking for
% plot c and p to see if we make sense
% investigate chebfun
% normalisation of chebyshev coefficients

for i = 1:m
    v = ((rand(1,dataset_size(1))<.5)*2 - 1)'; % randmacher vector;
    u = c(1)*v;
    if n>1 
        w0 = v;
        w1 = A_new*v;
        u = u + c(2)*w1;
        for j = 3:(n)
            w2 = 2*A_new*w1 - w0;
            u = u + c(j)*w2;
            w0 = w1;
            w1 = w2;
        end
    end
    fprintf('%.f\n', 100*i/m);
    G = G + v'*u/m;
end
% logdet(A_new)=logdet(alpha*A_old)
%              =log(det(alpha*I)*det(A_old))
%              =size*log(alpha)+logdet(A_old)
% logdet(A_old)=logdet(A_new)-size*log(alpha)
logdet = (G-dataset_size(1)*log(alpha))/2;
end