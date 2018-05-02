A_size = 8e4;
diag_dom_const = 1e-3;
while A_size <= 9e5
    epr = 2*floor(log10(A_size));
    A = Sparse_Dataset(A_size, epr, diag_dom_const);
    
    tic;
    amd = symamd(A);
    logdet_chol = full(2*sum(log(diag(chol(A(amd,amd))))));
    fprintf('it took %.3f seconds to find the logdet of a matrix of size %.f, and gave a value of %.4f\n', toc, A_size, logdet_chol);
    A_size = A_size+1e3;    
end