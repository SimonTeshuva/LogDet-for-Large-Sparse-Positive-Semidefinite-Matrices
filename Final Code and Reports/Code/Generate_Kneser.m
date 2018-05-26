%% Generate (n,k) Kneser Graph
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% this function takes 2 integers n and k, where n>=k and creates the (n,k)
% kneser graph. the edge list is then stored in a file 'kneser (n,k).txt'
% for later use.

% although this function can be used to generate any valid (n,k) Kneser
% Graph, only graphs of the form (2*k+1,k) were generated for this project.
% this is because those graphs guarantee the desired sparisity conitions
% for the experiments being run
%% Generate a Kneser Graph
function kneser_graph = Generate_Kneser(n,k)
% disjoint sets: no elements in common
% # vertices = n choose k
% # edges = (n choose k) * (n-k choose k)/2
% vertices correspond to k element subset of an n element set
% k element subset of an n element set: every combination of k elements in
% a set of size n
% each vertex has (n-k) choose k edges
no_vert = nchoosek(n,k);
no_edg = nchoosek(n,k)*nchoosek(n-k,k)/2;

degree = k+1;
% generate vertices
vertex_combos = nchoosek(1:n,k); % generates each k element subset of n.
% each subset represents a vertex
vertices = [(1:no_vert); vertex_combos']'; % assign each vertex a number 
% from 1 to n

% begin generating the edge list. compare each vertex with every other
% vertex. if they are disjoint (their k element subsets have no
% interseciton) they are adjacent, and an edge is added to the edge list
edges = [];
for i = 1:no_vert % iterating over every vertex pair in the graph to see if there is an edge
    vert1 = vertices(i, :);
    for j = i:no_vert
        vert2 = vertices(j, :);
        disjoint = true;
        for x = 2:length(vert1) % comparing the k element subsets of n of the two 
            % vertices. if disjoint, form an edge
            for y = x:length(vert2)
                if (vert1(x) == vert2(y)) || (vert1(y) == vert2(x))
                    disjoint = false;
                end
            end
        end
        if disjoint
            edges = [edges; [i, j]];
        end
    end
    if mod(i,10) == 0
        fprintf('%.2f\n', 100*i/no_vert);
    end
end

fprintf('this kneser graph has %.f vertices\n', no_vert);
fprintf('each of which has degree %.f\n', degree);
fprintf('resulting in %.knesef edges\n', no_edg);
fprintf('and an adjacency matrix of size %.f by %.f\n', no_vert, no_vert);
fprintf('with %.f elements in the full matrix\n', no_vert^2);

% for each edge in the edge list, print to a file for later use
fileName = strcat('kneser (',num2str(n),',',num2str(k),').txt');
fileID = fopen(fileName,'w');
for i = 1:length(edges)
    fprintf(fileID, '\n%.f %.f',edges(i,1), edges(i,2)); 
end

% also return the edge list as a sparse matrix for immediate use
kneser_graph = sparse(edges(:, 1), edges(:,2), 1, no_vert, no_vert);
end