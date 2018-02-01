function retVal = Experiment1
    density = 10;
    n = 15;
    m = 10;
    sub_size = 5;
    
    sigma_min = 1e-3;
    sigma_max = 1; % change to 1 norm of covariance;
    
    d_vect = [1e3, 3e3, 6e3, 1e4, 3e4, 6e4, 1e5, 3e5, 6e5, 1e6, 3e6, 6e6, 1e7, 3e7];
    d_vect_chol_schur = d_vect(1:sub_size);
%     d_vect = d_vect_chol_schur; % just for testing purposes
    
    logdet_vect = zeros(1, length(d_vect));
    logdet_vect_chol = zeros(1, sub_size);
    logdet_vect_schur = zeros(1, sub_size);
    
    t_vect = zeros(1, length(d_vect));
    t_vect_chol = zeros(1, sub_size);
    t_vect_schur = zeros(1, sub_size);
    
    % use Approx_Algorithm, Cholesky Decomposition, and Schur Complement
    
    % a) runtime of Approx_Algotithm for 1e3 --> 3e7
    
    for counter = 1:length(d_vect)
    fprintf('working on dataset of size %e for runtime of approx\n', d_vect(counter));
    
    dataset = Sparse_Dataset(d_vect(counter), density);
    
    tic
    logdet_vect(counter) = Approx_Algorithm(dataset, m, n);
    t_vect(counter) = toc;
    
    end
    
    figure(1)
    semilogy(t_vect, d_vect);
    title('runtime vs matrix dimension');
    xlabel('runtime');
    ylabel('matrix dimension');
    
    
    % b) relative accuracy for 1e3 --> 3e4    
    for counter = 1:length(d_vect_chol_schur)
        fprintf('working on dataset of size %e for relative accuracy\n', d_vect(counter));
        dataset = Sparse_Dataset(d_vect_chol_schur(counter), density);
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
    relative_accuracy = logdet_vect(1:sub_size)./logdet_vect_chol;
    figure(2)
    semilogy(relative_accuracy, d_vect_chol_schur);
	title('relative accuracy vs matrix dimension');
    xlabel('relative accuracy');
    ylabel('matrix dimension');

	% c) runtime for Cholesky, Schur, Approx_algorithm for 1e3 --> 3e4
    figure(3)
    semilogy(t_vect(1:sub_size), d_vect_chol_schur, t_vect_chol, d_vect_chol_schur, t_vect_schur, d_vect_chol_schur);
	title('runtime comparison for cholesky, schur and algorithm');
    xlabel('runtime');
    ylabel('matrix dimension');
    legend('algorithm', 'cholesky', 'schur');
    
    % d) compare accuracy to Zhang & Leithead, 2007. using n= 1000
    figure(4)

    aug_size = size(d_vect)-size(d_vect_chol_schur);
    aug_vect = zeros(1,aug_size(2));
    
    retVal = [d_vect; t_vect; [t_vect_chol aug_vect]; [t_vect_schur aug_vect]; [logdet_vect aug_vect]; [logdet_vect_chol aug_vect]; [logdet_vect_schur aug_vect]];
end