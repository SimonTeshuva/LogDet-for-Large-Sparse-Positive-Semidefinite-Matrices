%% Parameter Sweep for Various Variables
% Author: Simon Tesuva
% Date Last Modified: 14/5/2018

% This function uses the adjacency matrix of a kneser graph to sweep over 
% the condition number, function order, number of samples, estimates of
% lmin and lmax, and the number of times to run a method to average the
% results. the logdet is computed using the Chebyshev Function
% Appoximation, Rational Function Approximation, and Exact logdet. the
% results, as well as the compute times and the parameters used for each
% comprison are then written to an output file. 
function kappa_M_N_lmin_lmax(fileName, kappas, orders, samples,lambda_lows,lambda_highs,averages)
    %% Taking values for parameters
    if nargin == 0
        fileName = 'kneser (23,11).txt';
        kappas = [1.01, 1e10];
        orders = 14;
        samples = 20;
        lambda_lows = [1, 0.33, 0.01];
        lambda_highs = [1.75, 25, 1e3]; 
        averages = 3;    
    elseif nargin == 7
        fprintf('using values inputted form main');
    else
        return
    end
    
    %% generating the matrix and the file for results
    [n,k,A] = generate_adjacency_matrix(fileName);

    fileName = strcat('Ex3results(',num2str(n),',',num2str(k),').txt');
    fileID = fopen(fileName, 'a');
    
    parfor i=1:2
        fprintf('starting parallel pool\n');
    end
    
    
    %% run compare_methods on every permuation of the 6 parameters
    % and write results to output
    fprintf(fileID, 'results for experiment 3 at timestamp %s\n', strcat(num2str(fix(clock))));
    fprintf(fileID, 'samp#\torder\tkappa\taverages\texact\trational\ttime\tchebyshev\ttime\tlambda_min\tlambda_max\n');
    % for every combination of the inputs, run
    % compare_methods, which evaluates the logdet
    % using all three algorithms.
    for i = 1:length(samples)
        sample = samples(i);
        for j = 1:length(orders)
            order = orders(j);
            for x = 1:length(kappas)
                kappa = kappas(x);
                for y = 1:length(lambda_lows)
                    lambda_low = lambda_lows(y);
                    for z = 1:length(lambda_highs)
                        lambda_high = lambda_highs(z);
                        data_sum = 0;
                        for w = 1:averages % perform each compare_methods a number of times
                            % and average the results to increase
                            % reliablility of results
                            data_iteration = compare_methods(A,n,k,kappa,order,sample,lambda_low,lambda_high);
                            data_sum = data_sum+data_iteration;
                        end
                        data = data_sum/averages;
                        fprintf(fileID, '%.f\t\t%.f\t\t%.f\t\t%.f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.5f\t\t%.3f\t\t%.3f\n', sample, order, kappa, averages, data, lambda_low, lambda_high);
                    end
                end
            end
        end
        fprintf('we are %.2f percent of the way done\n',100*(i)/length(samples))
    end
    
    fclose(fileID);
end

%% Eigenvalues and Multiplicities for Kneser (n,k)

% compute the eigenvalues and their multiplicities eigenvectors 
% using known formulas
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

%% compare the three algorithms

% for a given matrix, and set of parameters being swept over, calculate the
% logdet using the exact formula, the Rational Function Method, and the
% Chebyshev Function Method, and store the result and the time taken to
% compute it.
function retVal = compare_methods(A,n,k,kappa,order,sampling_number,lmin,lmax)
    [evals, mults] = kneser_eigenvalues(n,k);
    I = speye(sum(mults), sum(mults));
    epsilon = (max(evals)-min(evals))/(kappa-1);
    
    logdet_exact = sum(mults.*log(evals - min(evals) + epsilon));
    
    A_new = diag(A) + (epsilon - min(evals))*I;
    
    tic
    logdet_rational = Rational(A_new,order,sampling_number,sampling_number);
    logdet_rational_time = toc;
    
    pertubated_lmin = lmin*epsilon;
    pertubated_lmax = lmax*(max(evals)-min(evals)+epsilon);
    
    tic
    logdet_chebyshev = ChebLogDet(A_new,sampling_number,order,pertubated_lmin,pertubated_lmax);
    logdet_chebyshev_time = toc;
    
    retVal = [logdet_exact, logdet_rational, logdet_rational_time, logdet_chebyshev, logdet_chebyshev_time];
end