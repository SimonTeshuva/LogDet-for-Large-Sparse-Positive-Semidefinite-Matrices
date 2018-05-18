% If you label the variables in your m (rows) x n (columns) grid down columns first (in other words, the correspondence is the labelling given by running the code reshape(1:m*n,m,n)) then 
% adj_path = diag(ones(m-1,1),1) + diag(ones(m-1,1),-1);
% will make the adjacency matrix of a length m path and
% adj_cycle = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1) + diag(1,n-1) + diag(1,-n+1);
% will make the adjacency matrix of a cycle with n vertices and
% 
% kron(adj_cycle,eye(m)) + kron(eye(n),adj_path)
% 
% will give the adjacency matrix of a graph that is an mxn grid with the (i,1) entry connected to the (i,n) entry for i=1,2,..,m
% 
% Please verify this on a small example before proceeding!

simple = 0;
if simple == 1
    m=5;n=5;
    adj_path=diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
    adj_cycle=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1);
    A = kron(adj_cycle,speye(m))+kron(speye(n),adj_path);
    spy(A)
else
    data=load('NASA_Data_needed.txt');
    m=length(data(:,2));
    n=length(data(:,3));

%     spdiags(data,d,m,n)
    adj_path=diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
    adj_cycle=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1);
    A = kron(adj_cycle,speye(m))+kron(speye(n),adj_path);
    spy(A)
    
end
