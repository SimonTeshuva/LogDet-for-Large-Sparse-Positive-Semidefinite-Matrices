function [A,long,lat,vals] = NASA_Exp(cbf)
    if nargin==0
        cbf = 0;
    elseif nargin == 1
        fprintf('using user input')
    else
        fprintf('invalid input');
        return;
    end

    if cbf
        [long,lat,vals] = read_data();
        m=length(long);n=length(lat);
        A = generate_adj(m,n);
        spy(A);
    else
        m=5;n=5;
        A = generate_adj(m,n);
        [long,lat,vals] = read_data();    
        for i = 1:length(long)
            row = round(100*long(i));
            col = round(100*lat(i));
            B(row,col) =  vals(i);
            fprintf('%.f%%\n',100*i/length(long));
        end
    end
end

function A = generate_adj(m,n)
    % reshape(1:m*n,m,n);
    adj_path = diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
    adj_cycle = diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1);

    A = sparse(kron(adj_cycle,speye(m)))+sparse(kron(speye(n),adj_path));
end

function [longditude,latitude,values] = read_data()
    data = load('NASA_Data.txt');
    longditude = data(:,2);
    latitude = data(:,3);
    values = data(:,6);
end