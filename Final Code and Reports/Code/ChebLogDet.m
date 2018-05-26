%% Chebyshev Function Appoximation for the Logdet

% Authors: I. Han, D. Malioutov, J. Shin. 
% “Large-scale Log-determinant Computation through Stochastic Chebyshev Expansions”, 

function ld = ChebLogDet( C, M, n, lmin, lmax )
% C : SPSD matrix with eigenvalues in [lmin , lmax]
% M : number of sampling for trace estimator
% n : Chebyshev approximate degree

if nargin == 4
    a = 1;
    delta = lmin;    
elseif nargin == 5
    a = lmin + lmax;
    delta = lmin / a;
end

% normalize B with eigenvalues in [delta, 1-delta]
B = C./a;

f = @(x) log(1-x);
g = @(x) ((1-2*delta)/2).*x+0.5;
ginv = @(x) (2/(1-2*delta)).*x;
h = @(x) f(g(x));

c = chebpolfit(h,n);

v = sign(randn(size(B,1),M));

u = c(1)*v;
if n>0
    w0 = v;
    w1 = (B*v);
    w1 = ginv(w1);
    w1 = v./(1-2*delta)-w1;
    u = c(2)*w1 + c(1)*w0;

    for j = 3 : n+1
        ww = (B*w1);
        ww = ginv(ww);
        ww = w1./(1-2*delta)-ww;
        ww = 2*(ww) - w0;
        u = c(j)*ww + u;
        w0 = w1;
        w1 = ww;
    end
end
ld = sum(sum(v.*u))/M + size(C,1) * log(a);
end

function c = chebpolfit(fname,n)
x = cos(((0:n)'+0.5)*pi/(n+1));
y = feval(fname,x);
T = [zeros(n+1,1) ones(n+1,1)];
c = [sum(y)/(n+1) zeros(1,n)];
a = 1;
for k = 2 : n+1
    T = [T(:,2) a*x.*T(:,2)-T(:,1)];
    c(k) = (y'*T(:,2))*2/(n+1);
    a = 2;
end
end