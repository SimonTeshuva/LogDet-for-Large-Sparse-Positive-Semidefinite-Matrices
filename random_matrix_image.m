A=Sparse_Dataset(1e2,5,1e-3);
colormap('gray');

figure(1)
spy(A);
title('Randomly Generated Diagonally Dominant Positive Definite Matrix');

[n,k,A]=generate_adjacency_matrix('kneser (7,3).txt');
kappa = 5;
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

% formula for value to add to diagonal
epsilon = (max(evals)-min(evals))/(kappa - 1);
lmin=min(evals); lmax=max(evals);
A=diag(A)+speye(size(A))*(epsilon-lmin);

figure(2)
imagesc(A);
title('Adjacency Matrix of (7,3) Kneser Graph with Condition Number set to 5');
