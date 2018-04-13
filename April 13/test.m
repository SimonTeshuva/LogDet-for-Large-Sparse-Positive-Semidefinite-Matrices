function kneser_graph = test(n,k)
% disjoint sets: no elements in common
% # vertices = n choose k
% # edges = (n choose k) * (n-k choose k)/2
% vertices correspond to k element subset of an n element set
% k element subset of an n element set: every combination of k elements in
% a set of size n
% each vertex has (n-k) choose k edges


no_vert = nchoosek(n,k);
no_edg = nchoosek(n,k)*nchoosek(n-k,k)/2;
degree = no_edg/no_vert;
% generate vertices
vertex_combos = nchoosek(1:n,k);
vertices = [(1:no_vert); vertex_combos']';
edges = [];
for i = 1:no_vert
    vert1 = vertices(i, :);
    for j = 1:no_vert
        vert2 = vertices(j, :);
        if sum(vert1(2:end) == vert2(2:end)) == 0
            edges = [edges; [i,j]];
        end
    end
end

fprintf('this kneser graph has %.f vertices\n', no_vert);
fprintf('each of which has degree %.f\n', degree);
fprintf('resulting in %.f edges\n', no_edg);
fprintf('and an adjacency matrix of size %.f by %.f\n', no_vert, no_vert);


kneser_graph = edges;
end