%% Generate Adjacency for (n,k) Kneser Graph
% Author: Simon Tesuva
% Date Last Modified: 8/4/2018

% this function takes the name of a file which stores a kneser graph and
% returns the adjacency matrix that describe it, as well as the kneser
% graph's n and k values

%% Generate the Adjacency Matrix
function [n,k,adjacency_matrix] = generate_adjacency_matrix(fileName)
    fileID = fopen(fileName,'r');
    tline = fgetl(fileID);
	% get the n and k values for the graph from the file name.
    % since all files are named in the format: 'kneser (n,k).txt' the
    % first 8 and last 5 characters can be ignored. the n and k can then be
    % found by spliting the remaining characters where the ',' is. 
    nk = split(fileName(9:(end-5)),',');    
    n = str2num(nk{1});
    k = str2num(nk{2});
    % kneser graphs have n choose k vertices
    no_vert = nchoosek(n,k);
    no_edge = nchoosek(n,k)*nchoosek(n-k,k)/2;
    edge_list = ones(no_edge,2);
    % read the edges in one at a time and add them to a growing edge list
    counter = 1;
    while ischar(tline)     
            tline=fgetl(fileID);
        if tline ~= -1
            edge = str2num(tline);
            edge_list(counter,:) = [edge(1), edge(2)];
            counter = counter+1;
        end
        if mod(counter,1e4)==0
            fprintf('%.3f\n',100*counter/no_edge)
        end
    end
	edge_list(counter,:) = [edge(1), edge(2)];

    % create a sparse matrix of size nchoosek x nchoose k and fill each
    % entry which has an edge with a 1
    A = sparse(edge_list(:,1)',edge_list(:,2)',ones(length(edge_list(:,1)),1),no_vert,no_vert);
    fclose(fileID);
    % the kneser matrix is undirected so the adjacency matrix needs to be
    % symetric. 
    adjacency_matrix=A|A';
end