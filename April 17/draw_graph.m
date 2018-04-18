function draw_graph(fileName)
    [n,k,A]=generate_adjacency_matrix(fileName);
    g=graph(A);
    titlestr=strcat('The (',num2str(n),',',num2str(k),') Kneser Graph');

    fig=figure;
    plot(g);
    title(titlestr);
    print(fig,titlestr,'-dpng');
end