function [logdet, results] = Parallel_Rational_Reordered(dataset, M, N, k)
    B = dataset;
    dataset_size = size(B);
    n = dataset_size(1);
    I = speye(n,n);
    [t, w] = lgwt(M, 0, 1);
	V = ((rand(N+k,n)<.5)*2 - 1)'; % randmacher vector;
    
    % logdet=(1/N)*SUM(1:N){v(i)'*(B-I)*SUM(1:M){wj*(tj*B+(1-tj)*I)^(-1)*v(i)}}
    % logdet
    % =SUM(1:M){(wj/N)*SUM(1:N){v(i)'*(B-I)*pcg(tj*B+(1-tj)*I),v(i)}}
    % =SUM(1:M){(wj/N)*SUM(1:N){v(i:(i+k))'*(B-I)*pcg(tj*B+(1-tj)*I),v(i:(i+k))}}

    gtotal = 0;
    parfor j = 1:M
        tj = t(j);
        wj = w(j);
        total = zeros(N+k, k);
        for i = 1:N
            v = V(:, i:(i+k-1));
            A = tj*B - (1-tj)*I;
            afun = @(v) reshape(A*(reshape(v, n, k)),n*k, 1);
            pcg_ret = reshape(pcg_quiet(afun, v(:)), n, k);
            LS_solved = wj*v'*(B-I)*pcg_ret;
            total(i:(i+k-1), :) = total(i:(i+k-1), :) + LS_solved;
        end
        adj_total = total(1:N, :);
        adj_total(1:k, :) = adj_total(1:k, :) + total((N+1):(N+k), :);        
        gtotal = gtotal+adj_total;
        fprintf('%.f\n', 100*j/M);
    end     
    results = gtotal;
    logdet = sum(sum(gtotal))/(N*k); 
end