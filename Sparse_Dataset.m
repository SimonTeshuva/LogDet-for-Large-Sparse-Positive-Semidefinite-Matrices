function retVal = Sparse_Dataset(size, epr, diag_dom_const)


% if we wand n elements per row, and the matrix is symetric, we need to
% generate n/2 points per row
no_samples = floor(epr/2);
% sparse matrix is stored as a list of ((i,j)    val)
fprintf('start sparse_dataset\n');
% ensure each collumn has the right number of elements by generating a
% vector in the form [1,1,1,1...2,2,2,2....3,3,3,3....] using kron
i = kron((1:size)', ones(no_samples, 1));
fprintf('i\n');
% assign each i point in each column with a random j point
j = randi(size, [no_samples*size, 1]);
fprintf('j\n');
% assign each i,j pair a random number in the range [-1,1]
val = -1+2*rand([no_samples*size, 1]);
fprintf('val\n');
% store as a sparse matrix
S = sparse(i, j, val);
fprintf('S_1\n');
% to make the matrix symetric, S+S'. subtract diagonal of S to take into
% account diagonal elements being doubled by this process
S = S+S'-2*diag(diag(S));
fprintf('S_2\n');
% make the matrix diagonally dominant by adding setting the diagonals to be
% 1e3+abs(sum(row))
val_d = sum(abs(S'))+diag_dom_constant;
% S_diag = sparse(1:size, 1:size, val_d);
% S = S + S_diag;
S = spdiags(val_d', 0, S);

fprintf('S_DD\n');

retVal = S;

end

% 1e4 is instant and is 2MB large
% 1e5 takes a second and is 20MB large
% 1e6 takes 5 seconds and is 200MB large
% 1e7 takes 1 minute and is arouind 2GB large