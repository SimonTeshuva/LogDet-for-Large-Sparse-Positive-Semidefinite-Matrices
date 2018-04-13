function logdet = Rational(dataset, M, N, k)
% optimal N is 20 done in 1 block

B = dataset;
dataset_size = size(B);
n = dataset_size(1);
I = speye(n,n);
[t, w] = lgwt(M, 0, 1);
V = ((rand(N,n)<.5)*2 - 1)'; % randmacher vector;

% logdet=(1/N)*SUM(1:N){v(i)'*(B-I)*SUM(1:M){wj*(tj*B+(1-tj)*I)^(-1)*v(i)}}
% logdet
% =SUM(1:M){(wj/N)*SUM(1:N){v(i)'*(B-I)*pcg(tj*B+(1-tj)*I),v(i)}}
% =SUM(1:M){(wj/N)*SUM(1:N){v(i:(i+k))'*(B-I)*pcg(tj*B+(1-tj)*I),v(i:(i+k))}}
fake_N = k*floor(N/k);
gtotal = 0;
parfor j = 1:M
    tj = t(j);
    wj = w(j);
    total = 0;
    for i = 1:k:fake_N
        v = V(:, i:(i+k-1));
        A = tj*B + (1-tj)*I;
        afun = @(v) reshape(A*(reshape(v, n, k)),n*k, 1);
        pcg_ret = reshape(pcg_quiet(afun, v(:)), n, k);
        BIpcg = (B-I)*pcg_ret;
        
        LS_solved = wj*v(:)'*BIpcg(:);
        total = total + LS_solved;
    end
    gtotal = gtotal + total;
    
    fprintf('%.f\n', 100*j/M);
end
logdet = gtotal/fake_N;
%     logdet = sum(sum(gtotal))/(N*k);
end