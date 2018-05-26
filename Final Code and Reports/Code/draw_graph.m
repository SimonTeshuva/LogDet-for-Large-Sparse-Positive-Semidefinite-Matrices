%% Draw Kneser Graph 
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% this function takes as input the name of a file containing the edge list
% of an (n,k) kneser graph. It then reads from the file name what the n and
% k values are, plots the Kneser Graph, and prints it to a file titled 
% 'The (n,k) Kneser Graph.png', where n and k are the values found from the
% file name
function draw_graph(fileName)
    [n,k,A]=generate_adjacency_matrix(fileName);
    g=graph(A);
    titlestr=strcat('The (',num2str(n),',',num2str(k),') Kneser Graph');

    fig=figure;
    plot(g);
    title(titlestr);
    print(fig,titlestr,'-dpng');
end