%% Generate Reduced Laplacian
% Author: Simon Tesuva
% Date Last Modified: 8/4/2018

% this function takes in the name for a file coressponding with a (n,k)
% kneser graph returns a sparse matrix coresponding to the reduced
% laplacian RL(A) of that matrix.

% the full laplacian L(A) is first found:
% L(A) = D-A
% where D = a diagonal matrix corresponding to the degree of each vertex.
% in the case of a (2*k+1, k) kneser graph, the graph is regular and each
% vertex has degree k+1, so L(A) can be expressed as;
% L(A) = (k+1)*I-A
% otherwise, D must be found by taking the sum of each row to find the
% degree of each vertex. 
% in the case of the (2*k+1,k) graph, both formulas will return the same
% result, but the first formula will involve far fewer calculations.

% the reduced laplacian is found by ignoring the row and column
% corresponding to the only eigenvalue multiplicity of 1, effectivly
% removing its vertex from the graph. 
% Since the kneser graph is regular, any vertex can be removed. for
% simplicity, this function always removes the last vertex.
%% Code
function reduced_laplacian = Laplacian(fileName)
    [n,k,A] = generate_adjacency_matrix(fileName);
    if n == 2*k+1
        D = (k+1)*speye(size(A));
    else
        D = diag(sum(A));
    end
    laplacian = D-A;
    reduced_laplacian=laplacian(1:(end-1),1:(end-1));
end