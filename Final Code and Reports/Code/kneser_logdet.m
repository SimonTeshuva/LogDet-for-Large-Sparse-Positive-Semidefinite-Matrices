%% Logdet of Condition Number Adjusted Kneser Graph Adjacency Matrix
% Author: Simon Tesuva
% Date Last Modified: 13/5/2018

% This function takes input the n and k values for a kneser graph, as well
% as a condition number kappa, and returns the logdet of the resultant
% matix, as well as the condition number shifted eigenvalues and their
% multiplicities
function [logdet, eigenvalues, multiplicities] = kneser_logdet(n,k,kappa)
    evals = [];
    mults = [];
    for i = 0:k
        % formula for the eigenvalues of a kneser graph matrix
        eval = ((-1)^i)*nchoosek(n-k-i, k-i); 
        % formula for the multiplicities of each eigenvalue of a kneser
        % graph matrix
        if i == 0
            mult = 1;
        else 
            mult = nchoosek(n,i) - nchoosek(n, i-1);
        end
        evals = [evals eval];
        mults = [mults mult];
    end
    
    % formula for value to add to diagonal to set the condition number to
    % kappa
	epsilon = (max(evals)-min(evals))/(kappa - 1);
    evals = evals - min(evals) + epsilon;
    
    % logdet is equivaluent to taking the sum of the log of the eigenvalues
    logdet = sum(mults.*log(evals));
    eigenvalues = evals;
    multiplicities = mults;
end