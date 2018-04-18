function retVal = Experiment3(fileName, kappa_min, kappa_max, kappa_checknum, function_order_min, function_order_max, sampling_number_min, sampling_number_max, sampling_number_checknum)
    [n,k,A] = generate_adjacency_matrix(fileName);
    
    fileID = fopen('Ex3results.txt', 'a');
    
    parfor kappa_counter=1:2
        fprintf('starting parallel pool\n');
    end
    
    fprintf(fileID, 'results for experiment 3 at timestamp %s\n', strcat(num2str(fix(clock))));
    fprintf(fileID, 'samp#\t\torder\t\tkappa\t\texact\t\trational\t\ttime\t\tchebyshev\t\ttime\n');
    
    for sampling_number_counter = 0:(sampling_number_checknum)
        sampling_number = sampling_number_min + sampling_number_counter*((sampling_number_max-sampling_number_min)/(sampling_number_checknum));
        for order = function_order_min:function_order_max
            for kappa_counter = 0:kappa_checknum
                kappa = kappa_min + kappa_counter*(kappa_max - kappa_min)/(kappa_checknum);
                data = compare_methods(A,n,k,kappa,order,sampling_number);
                fprintf(fileID, '%.0f\t\t%.0f\t\t%.2f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\n', sampling_number, order, kappa, data);
            end
        end
        fprintf('we are %.2f percent of the way done\n',100*(sampling_number_counter+1)/sampling_number_checknum)
    end
    
    fclose(fileID);
end

function [eigenvalues, multiplicities] = kneser_eigenvalues(n,k)
    evals = [];
    mults = [];
    for i = 0:k
        eval = ((-1)^i)*nchoosek(n-k-i, k-i); 
        if i == 0
            mult = 1;
        else 
            mult = nchoosek(n,i) - nchoosek(n, i-1);
        end
        evals = [evals eval];
        mults = [mults mult];
    end
    eigenvalues = evals;
    multiplicities = mults;
end

function retVal = compare_methods(A,n,k,kappa,order,sampling_number)
    [evals, mults] = kneser_eigenvalues(n,k);
    I = speye(sum(mults), sum(mults));
    epsilon = (max(evals)-min(evals))/(kappa-1);
    
    logdet_exact = sum(mults.*log(evals - min(evals) + epsilon));
    
    A_new = A + (epsilon - min(evals))*I;
    
    tic
    logdet_rational = Rational(A_new,order,sampling_number,sampling_number);
    logdet_rational_time = toc;
    
    tic
    logdet_chebyshev = ChebLogDet(A_new,sampling_number,order,epsilon,max(evals)-min(evals)+epsilon);
    logdet_chebyshev_time = toc;
    
    retVal = [logdet_exact, logdet_rational, logdet_rational_time, logdet_chebyshev, logdet_chebyshev_time];
end