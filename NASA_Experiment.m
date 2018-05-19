% function NASA_Experiment(make_data)
%     if nargin == 0
%         make_data = 0;
%     elseif nargin == 1
%         fprintf('using user input\n');
%     else
%         return
%     end
    make_data=0;
    if make_data == 1
        fileName = 'NASA_Data.txt';
        data_full = load(fileName);
        long = data_full(:,2); % long and lat may be the wrong way around;
        lat = data_full(:,3);
        val = data_full(:,6);
        data_needed = [long,lat,val];
        clear data_full long lat val;
        fileName = 'NASA_Data_needed.txt';
        fileID = fopen(fileName,'w+');
        for i=1:size(data_needed,1)
            fprintf(fileID, '%.2f %.2f %.2f\n',data_needed(i,:));
            if mod(i,5e3)==0
                fprintf('%%%.3f\n',100*i/size(data_needed,1));
            end
        end
        fclose(fileID);
    else
        data_needed = load('NASA_Data_needed.txt');
    end
    
    fprintf('starting\n');
    integer_vals_long = round(100*(data_needed(:,1)-min(data_needed(:,1))))+1;
    integer_vals_lat = round(100*(data_needed(:,2)-min(data_needed(:,2))))+1;
    Y=data_needed(:,3);
    
    fprintf('got data\n');
    
    save('temp.mat');
    % probably can be done better with find()
    m=max(max(integer_vals_long),180*100*2);
    n=max(max(integer_vals_lat),90*100*2);
    
    fprintf('got Mobs\n');
    Mobs=logical(sparse(integer_vals_long,integer_vals_lat,Y,m,n));
    
    z=reshape(1:m*n,m,n);
    M_obs_indecies = z(Mobs);
	M_hidden_indicies = z(~Mobs);
    fprintf('got Mobs indices\n');

    adj_path=sparse(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1));
    fprintf('got adj_path\n');
    adj_cycle=sparse(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1));
    fprintf('got adj_cycle\n');
    
	right = Mobs(:,[2:end,1]);
    left = Mobs(:,[end,1:end-1]);
    down = [Mobs(2:end,:);zeros(1,n)];
    up = [zeros(1,n);Mobs(1:end-1,:)];
    deg_mat=right+left+down+up;
    fprintf('got degrees of vertices\n');
    deg=reshape(deg_mat,m*n,1);
    
    cycle = kron(adj_cycle,speye(m));
    path = kron(speye(n),adj_path);
    A=cycle+path;
    fprintf('got A\n');
    
    Jtm = spdiag(deg)-A;
    V=Jtm*diag(1./deg);
    Jtp=V*V';
    
    alpha_vector = -10:2:10;
    beta_vector = -10:2:10;
    
    [alpha,beta]=max_alpha_beta(alpha_vect,beta_vect,Y,Jtp,Jtm);
    
    I=speye(size(Jtm));
    
	J = alpha*I + beta*Jtp + (1-beta)*Jtm;
    
% 	Jxx = J(M_obs_indecies,M_obs_indecies);
%	Jxy = J(M_obs_indecies,M_hidden_indicies);
% 	Jyx = J(M_hidden_indicies,M_obs_indecies);
%	Jyy = J(M_hidden_indicies,M_hidden_indicies);
    
    X=pcg(J(M_obs_indecies,M_obs_indecies),J(M_obs_indecies,M_hidden_indicies)*Y);
    
    imshow(Y,X);

function [alpha_max, beta_max] = max_alpha_beta(alpha_vect,beta_vect,Y,Jtp,Jtm)
    alpha_beta_matrix = zeros(length(alpha_vector),length(beta_vector));
    I = speye(size(Jtp));
    for i=1:length(alpha_vect)
        alpha=alpha_vect(i);
        for j=1:length(beta_vect)
            beta=beta_vect(j);
            J = alpha*I + beta*Jtp + (1-beta)*Jtm;
    
            logdet_J = Rational(J,14,20,20); % if this takes too long, down sample
            logdet_Jxy = Rational(J(M_obs_indecies,M_hidden_indicies),14,20,20);

            loglike_J=0.5*logdet_J;
            loglike_J=loglike_J-0.5*logdet_Jxy;
            loglike_J=loglike_J-0.5*Y'*J(M_hidden_indicies,M_hidden_indicies)*Y;
            loglike_J=loglike_J+0.5*Y'*J(M_hidden_indicies,M_obs_indecies)*pcg(J(M_obs_indecies,M_obs_indecies),J(M_obs_indecies,M_hidden_indicies)*Y);
            
            alpha_beta_matrix(i,j)=loglike_J;
        end
    end
    
    maximized_loglike = max(max(loglike_J));
    alpha_beta = index(maximized_loglike); 
    alpha_max=alpha_beta(1);
    beta_max=alpha_beta(2);
end