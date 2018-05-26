% take as input the k value for a Kneser Graph and return the file name of
% the associated (n,k) Kneser Graph, where n=2*k+1
function fileName = kneser_file(k)
    fileName=strcat('kneser (',num2str(2*k+1),',',num2str(k),').txt');
end