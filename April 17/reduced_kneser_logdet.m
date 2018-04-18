function [logdet, mults, evals] = reduced_kneser_logdet(n,k)
% http://math.mit.edu/~csikvari/spectral_graph_theory_V3.5.pdf
    eigenvalues = zeros(1,k);
    multiplicities = zeros(1,k);
    for i = 0:k
        eigenvalues(i+1)=(-1)^(i)*(nchoosek(n-k-i,k-i));
        if i==0
            multiplicities(i+1)=1;
        else
            multiplicities(i+1)=nchoosek(n,i)-nchoosek(n,i-1);
        end
    end
    mltlev = multiplicities.*log(eigenvalues);
    logdet=2*sum(mltlev(2:end));
    mults=multiplicities;
    evals=eigenvalues;
end