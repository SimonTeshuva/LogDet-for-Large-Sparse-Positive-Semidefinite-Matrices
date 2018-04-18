function [logdet, eigenvalues, multiplicities] = kneser_logdet(n,k,kappa)
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

	epsilon = (max(evals)-min(evals))/(kappa - 1);
    evals = evals - min(evals) + epsilon;
    
    logdet = sum(mults.*log(evals));
    eigenvalues = evals;
    multiplicities = mults;
end