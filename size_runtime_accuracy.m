function retVal = size_runtime_accuracy
    epr = 10;
    diag_dom_const = 1e-3;
    n = 15;
    m = 10;
    sub_size = 5;
    k = min(20, m);
    

    d_vect = [1e3, 3e3, 6e3, 1e4, 3e4, 6e4, 1e5, 3e5, 6e5, 1e6, 3e6, 6e6];
        %, 1e7, 3e7];
    d_vect_chol_schur = d_vect(1:sub_size);
%     d_vect = d_vect_chol_schur; % just for testing purposes
    
    logdet_vect_rat = zeros(1, length(d_vect));
    logdet_vect_cheb = zeros(1, length(d_vect));
    logdet_vect_chol = zeros(1, sub_size);
    logdet_vect_schur = zeros(1, sub_size);

	t_vect_rat = zeros(1, length(d_vect));
    t_vect_cheb = zeros(1, length(d_vect));
    t_vect_chol = zeros(1, sub_size);
    t_vect_schur = zeros(1, sub_size);
    
    % use Approx_Algorithm, Cholesky Decomposition, and Schur Complement
    
    % a) runtime of Approx_Algotithm for 1e3 --> 3e7
    parfor i=1:2
        fprintf('Experiment 1');
    end
    
    for counter = 1:length(d_vect)
        fprintf('working on dataset of size %e for runtime of approx\n', d_vect(counter));

        dataset = Sparse_Dataset(d_vect(counter), epr, diag_dom_const);
        
        sigma_min = 1e-3;
        sigma_max = max(max(dataset)); % change to 1 norm of covariance;

        tic
        logdet_vect_cheb(counter) = ChebLogDet(dataset, m, n, sigma_min, sigma_max);
        t_vect_cheb(counter) = toc;
        
        tic
        logdet_vect_rat(counter) = Rational(dataset, n, m, k);
        t_vect_rat(counter) = toc;        
    end
    
    figure(1)
    semilogy(t_vect_cheb, d_vect, t_vect_rat, d_vect);
    title('runtime vs matrix dimension');
    xlabel('runtime');
    ylabel('matrix dimension');
    legend('chebyshev', 'rational');
    
    % b) relative accuracy for 1e3 --> 3e4    
    for counter = 1:length(d_vect_chol_schur)
        fprintf('working on dataset of size %e for relative accuracy\n', d_vect(counter));
        dataset = Sparse_Dataset(d_vect_chol_schur(counter), epr, diag_dom_const);
        tic
        R = chol(dataset);
        logdet_vect_chol(counter) = 2*sum(log(diag(R)));
        t_vect_chol(counter) = toc; 
        fprintf('cholesky done\n');
        
        tic
        T = schur(full(dataset));
        logdet_vect_schur(counter) = 2*sum(log(diag(T))); % not sure if this is what i should be doing;
        t_vect_schur(counter) = toc;
        fprintf('schur done\n');
        
    end
    relative_accuracy_cheb = (logdet_vect_chol - logdet_vect_cheb(1:sub_size))./logdet_vect_chol;
    relative_accuracy_rat = (logdet_vect_chol - logdet_vect_rat(1:sub_size))./logdet_vect_chol;
    figure(2)
    semilogy(relative_accuracy_cheb, d_vect_chol_schur, relative_accuracy_rat, d_vect_chol_schur);
	title('relative accuracy vs matrix dimension');
    xlabel('relative accuracy');
    ylabel('matrix dimension');
    legend('chebyshev', 'rational');

	% c) runtime for Cholesky, Schur, Approx_algorithm for 1e3 --> 3e4
    figure(3)
    semilogy(t_vect_cheb(1:sub_size), d_vect_chol_schur, t_vect_chol, d_vect_chol_schur, t_vect_schur, d_vect_chol_schur, t_vect_rat(1:sub_size), d_vect_chol_schur);
	title('runtime comparison for cholesky, schur and algorithm');
    xlabel('runtime');
    ylabel('matrix dimension');
    legend('chebyshev', 'cholesky', 'schur', 'rational');
    
    % d) compare accuracy to Zhang & Leithead, 2007. using n= 1000
    figure(4)

    aug_size = size(d_vect)-size(d_vect_chol_schur);
    aug_vect = zeros(1,aug_size(2));
    
    retVal = [d_vect  ; t_vect_cheb; [t_vect_chol aug_vect]; [t_vect_schur aug_vect]; [logdet_vect_cheb aug_vect]; [logdet_vect_chol aug_vect]; [logdet_vect_schur aug_vect]];
end