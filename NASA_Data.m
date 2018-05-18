

function retVal = NASA_Data(ideal_sampling_number, ideal_function_order)
    A = generate_matrix_NASA();
    sizeA = size(A);
    n = sizeA(1);
    maxsize = 3e4;
    
    if n<maxsize
        A_reordered = greedy_reordering(A);
        tic; logdet_exact = 2*sum(log(diag(chol(A_reordered)))); logdet_exact_time = toc;
        fprintf('exact val is %.3f and took %.3f seconds to compute', logdet_exact, logdet_exact_time);
    end
    
    N = ideal_sampling_number;
    M = ideal_function_order;
    k = 20;
    eigs = get_eigenvalues(A);
    lambda_min = min(eigs);
    lambda_max = max(eig);
    tic; logdet_cheb = ChebLogDet(A,N,M,lambda_min,lambda_max); loget_cheb_time = toc;
	fprintf('exact val is %.3f and took %.3f seconds to compute', logdet_cheb, loget_cheb_time);

    tic; logdet_rational = Rational(A,M,N,k); logdet_rational_time = toc;
	fprintf('exact val is %.3f and took %.3f seconds to compute', logdet_rational , logdet_rational_time );    
    
    subplot(1,3,1), plot(A), title('unmodified data');
    subplot(1,3,2), plot(interpolate_data(A, logdet_cheb)), title('interpolated with chebyshev logdet');
    subplot(1,3,3), plot(interpolate_data(A, logdet_rational)), title('interpolated with rational logdet');
end

function A = generate_matrix_NASA()
    [n,k,A] = generate_adjacency_matrix('kneser (5,2).txt');
end

function A_reordered = greedy_reordering(A)
    A_reordered = A;
end

function retVal = interpolate_data(A, logdet)
    retVal = A - logdet+logdet;
end