%% NASA Data Interpolation Experiment
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% this experiment uses a dataset which describes the total worldwide ozone
% levels. the data is incomplete and the goal of this experiment is to show
% that the Rational Function Approximation Algorithm can be used as a key
% subroutine in performing a meaningful interpolation of the data. 

% the function takes as input a variable make_data, which determines if the
% dataset is generated from scratch or if it already exists. it should only
% be set to 1 the first time this experiment is run on any given machine.
function NASA_Experiment(make_data)

%% evaluating input
    if nargin == 0
        make_data = 0;
    elseif nargin == 1
        fprintf('using user input\n');
    else
        return
    end

%% obtaining data
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
 
    precision = .1; % data is rounded to this degree of precision to downsample it.
    % without downsampleing, the matrix that describes the dataset would be of size 
    % 6.5e8 x 6.5e8
    
    fprintf('starting\n');
    % map the logitude and latitude coordinates from doubles to integer values
    integer_vals_long = round((1/precision)*(data_needed(:,1)-min(data_needed(:,1))))+1;
    integer_vals_lat = round((1/precision)*(data_needed(:,2)-min(data_needed(:,2))))+1;
    Y=data_needed(:,3);
    
    fprintf('got data\n');

%% Generating Necessary Matrices

    m=180*2/precision+1; % longitude ranges in [-180,180] in steps of "1/precision"
    n=90*2/precision +1; % longitude ranges in [-180,180] in steps of "1/precision"
    
    % an m by n matrix. if there exists a datapoint at some coordinate,
    % place a 1 there. else, place a 0.
    Mobs=logical(sparse(integer_vals_long,integer_vals_lat,Y,m,n));
    fprintf('got Mobs\n');
    
    z=reshape(1:m*n,m,n);
    M_obs_indecies = z(Mobs); % a vector in column order of Mobs containing
    % all the observed datapoints
    fprintf('got Mobs indices\n');
    
	right = Mobs(:,[2:end,1]);
    left = Mobs(:,[end,1:end-1]);
    down = [Mobs(2:end,:);zeros(1,n)];
    up = [zeros(1,n);Mobs(1:end-1,:)];
    deg_mat=right+left+down+up;
    fprintf('got degrees of vertices\n');
    deg=reshape(deg_mat,m*n,1); % each coordinate is affected by its spacial 
    % neighbour to its left, right, up and down. shift the Mobs matrix by 1
    % in all these locations, and then add those matrices to get the degree
    % of each spacial co-ordinate. wrap around happens horizontally because
    % of how longitude wraps around [-180 = 180], but not around latitude as it does not
    % occur there [90 = north pole, -90 = south polt]
    
    adj_path=sparse(diag(ones(m-1,1),1)+diag(ones(m-1,1),-1)); % going down the latitude
    % lines there exist paths of vertices going from north to south
    fprintf('got adj_path\n');
    adj_cycle=sparse(diag(ones(n-1,1),1)+diag(ones(n-1,1),-1)+diag(1,n-1)+diag(1,-n+1));
    % going around the longitude there exist cycles of vertices wrapping
    % around the globe
    fprintf('got adj_cycle\n');
   
    cycle = kron(adj_cycle,speye(m)); % create a matrix that describes the cycles
    % using cron
    fprintf('got cycle\n');
    path = kron(speye(n),adj_path); % creat a matrix that describes the paths
    % using cron
    fprintf('got path\n');
    A=cycle+path; % the adjacency matrix is made up of the cycles and paths in the graph
    
    fprintf('got A\n');
    
    Jtm = sparse(diag(deg)-A); % thin membrane matrix used to approximate the information
    % matrix
	fprintf('got Jtm\n');

    V = sparse(Jtm.*(1./deg));
	fprintf('got V\n');
    Jtp=V*V'; % thin plate matrix used to describe the information matrix = 
    % V*V'
	fprintf('got Jtp\n');

    % course sweep over alpha and beta to start. then a fine sweep later to
    % get more accurate values
    alpha_vector = -10:2:10; 
    beta_vector = -10:2:10;
    
    % a vector in column order of Mobs containing all the un-observed datapoints
    M_hidden_indicies = z(~Mobs);
    fprintf('got M_hidden_indices\n');
    
    % maximize over alpha and beta max_alpha_beta
    [alpha,beta]=max_alpha_beta(alpha_vector,beta_vector,Y,Jtp,Jtm,M_obs_indecies,M_hidden_indicies);
    
    I=speye(size(Jtm));
    
    % generate the information matrix that describes the dataset
	J = alpha*I + beta*Jtp + (1-beta)*Jtm;
    
    
%   This is the logic being used to generate the Jxx,Jxy,Jyx & Jyy matrices
%   at the moment. it creates non-square matrices, so clearly something is
%   wrong. this is what needs to be debugged.
% 	Jxx = J(M_obs_indecies,M_obs_indecies);
% 	Jxy = J(M_obs_indecies,M_hidden_indicies);
% 	Jyx = J(M_hidden_indicies,M_obs_indecies);
% 	Jyy = J(M_hidden_indicies,M_hidden_indicies);

%   Jxx = the observed data points being affected by observed data points
%   Jxy = the observed data points being affected by un-observed data
%   points
%   Jyx = the un-observed data points being affected by observed data
%   points
%   Jyy = the un-observed data points being affected by un-observed data
%   points
	
%   interpolate the missing data using the formula shown in the report
    X=pcg(J(M_obs_indecies,M_obs_indecies),J(M_obs_indecies,M_hidden_indicies)*Y);
    
%   once the interpolation has been successfully done, this line will be swapped 
% with one that generates an image of the interpolated dataset.
    imshow(Y,X);

end

%% maximising alpha and beta
function [alpha_max, beta_max] = max_alpha_beta(alpha_vector,beta_vector,Y,Jtp,Jtm,M_obs_indecies,M_hidden_indicies)
    alpha_beta_matrix = zeros(length(alpha_vector),length(beta_vector));
    I = speye(size(Jtp));
    % sweep over alpha and beta and generate new matrices J to describe the
    % dataset, then evaluate the log-likelihood for each. the [alpha,beta] 
    % pair which maximises the log-likelihood is then returned to be used
    % to create the final J matrix which will be used in the interpolation
    for i=1:length(alpha_vector)
        alpha=alpha_vector(i);
        for j=1:length(beta_vector)
            beta=beta_vector(j);
            J = alpha*I + beta*Jtp + (1-beta)*Jtm;
            
% log-likelihood(J) =
% 0.5*logdet(J)-0.5*logdet(Jxy)-0.5*Y'*Jxx*Y+0.5*Y*Jyx*pcg(Jxx,Jxy*Y);
           

%             matrix dimension mismatch. clearly im indexing the matriceis
%             incorrectly. most steps from here onward.
            logdet_J = Rational(J,14,20,20);
            logdet_Jxy = Rational(J(M_obs_indecies,M_hidden_indicies),14,20,20);

            loglike_J=0.5*logdet_J;
            loglike_J=loglike_J-0.5*logdet_Jxy;
            loglike_J=loglike_J-0.5*Y'*J(M_hidden_indicies,M_hidden_indicies)*Y;
            loglike_J=loglike_J+0.5*Y'*J(M_hidden_indicies,M_obs_indecies)*pcg(J(M_obs_indecies,M_obs_indecies),J(M_obs_indecies,M_hidden_indicies)*Y);
            
            alpha_beta_matrix(i,j)=loglike_J;
        end
    end
    
    % find the alpha and beta values which maximised the log-likelihood
    % function and return them as output
    maximized_loglike = max(max(loglike_J));
    alpha_beta = index(maximized_loglike); 
    alpha_max=alpha_beta(1);
    beta_max=alpha_beta(2);
end