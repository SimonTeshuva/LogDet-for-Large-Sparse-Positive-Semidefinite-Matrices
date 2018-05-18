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
    
    integer_vals_long = round(100*(data_needed(:,1)-min(data_needed(:,1))))+1;
    integer_vals_lat = round(100*(data_needed(:,2)-min(data_needed(:,2))))+1;
    
    save('temp.mat');
    % probably can be done better with find()
    m=180*100*2;
    n=90*100*2;
    Mobs = logical(sparse(integer_vals_long,integer_vals_lat,1,m,n,length(integer_vals_long))); % last one not imporant
    
    z=reshape(1:m*n,m,n);
    M_obs_indecies = z(Mobs);
    M_hidden_indicies = z(~Mobs);
    
    Y=data_needed(:,3);
    
    Jtm = Generate_Jtm(data_needed);
    V = Generate_V(data_needed);
    Jtp = V*V';
    
    alpha_vector = -10:2:10;
    beta_vector = -10:2:10;
    
    [alpha_max, beta_max] = max_alpha_beta(alpha_vect,beta_vect,Y);
	I = speye(size(Jtp));

    J = alpha_max*I + beta_max*Jtp + (1-beta_max)*Jtm;
	Jxy = J;
    Jyy = J;
    
%     X=inv(Jxy)*(inv(inv(Jyy))*Y);
%     X=inv(Jxy)*(Jyy*Y); % need to not take inverse
    % Jxy*X=Jyy*Y
    X=pcg(Jxy,Jyy*Y); % needs updating
        
    % generate matrix B with Y & X in right locations
    % start with zeros() & fill with boolean index.
    imagesc(B) % may need to flip indexes if image is upside down/reversed (fliplr)
    
% end

function Jtm = Generate_Jtm(data)
%     adj_path=diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
%     adj_cycle=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1);
%     A = kron(adj_cycle,speye(m))+kron(speye(n),adj_path);
    Jtm = speye(length(data(:,1)),length(data(:,1)));
end

function V = Generate_V(data)
%     adj_path=diag(ones(m-1,1),1)+diag(ones(m-1,1),-1);
%     adj_cycle=diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1);
%     A = kron(adj_cycle,speye(m))+kron(speye(n),adj_path);
    V = speye(size(data));
end

function [alpha_max, beta_max] = max_alpha_beta(alpha_vect,beta_vect,Y)
    alpha_beta_matrix = zeros(length(alpha_vector),length(beta_vector));
    I = speye(size(Jtp));
    for i=1:length(alpha_vect)
        alpha=alpha_vect(i);
        for j=1:length(beta_vect)
            beta=beta_vect(j);
            J = alpha*I + beta*Jtp + (1-beta)*Jtm;
            Jxx = J;
            Jxy = J;
            Jyx = J;
            Jyy = J;
            logdet_J = Rational(J,14,20,20); % if this takes too long, down sample
            logdet_Jxy = Rational(Jxy,14,20,20);

            loglike_J=0.5*logdet_J;
            loglike_J=loglike_J-0.5*logdet_Jxy;
            loglike_J=loglike_J-0.5*Y'*Jyy*Y;
            loglike_J=loglike_J+0.5*Y'*Jyx*pcg(Jxx,Jxy*Y);
            
            alpha_beta_matrix(i,j)=loglike_J;
        end
    end
    
    maximized_loglike = max(max(loglike_J));
    alpha_beta = index(maximized_loglike); 
    alpha_max=alpha_beta(1);
    beta_max=alpha_beta(2);
end