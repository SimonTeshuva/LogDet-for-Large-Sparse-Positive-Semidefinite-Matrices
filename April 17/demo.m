function demo
N = 10000;
rho = 10/N;
%delta = 0.01;
%B = make_scaled_spd_matrix(N, rho, delta);

lmin = 0.1;

A = make_spd_matrix( N, rho, lmin);
lmax = norm(A,inf);
try
    tic;
    ld = 2*sum(log(diag(chol(A))))+0;
    fprintf('exact log-det : %f , elapsed time : %f\n', ld,toc);
catch err1
    fprintf('calculate exact log-det error.\n');
end

try
    tic;
    ld2 = ChebLogDet(A,15,40,lmin,lmax);
    fprintf('approx log-det : %f , elapsed time : %f\n', ld2,toc);
catch err2
    fprintf('calculate approx log-det error.\n');
end

fprintf('error rate : %f\n', abs(ld-ld2)/ld);
end

function A = make_spd_matrix(N, rho, lmin)
A = sprandsym(N,rho);
A = A * A';
A = A + lmin*speye(N);
end


function B = make_scaled_spd_matrix(N, rho, delta)

lambda_min = delta;
lambda_max = 1-lambda_min;

A = sprandsym(N,rho);

lmin = eigs(A,1,'SA');
lmax = eigs(A,1,'LA');

rho_corr = (lambda_max - lambda_min) / (lmax-lmin);
B = A*rho_corr + (lambda_min - lmin * rho_corr) *speye(N);
end

