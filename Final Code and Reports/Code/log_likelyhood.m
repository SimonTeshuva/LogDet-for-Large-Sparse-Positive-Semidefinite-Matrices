%% Log-Likelihood Maximisation
% Author: Simon Tesuva
% Date Last Modified: 26/5/2018

% this function takes as input the name of a kneser graph, a hidden value
% for rho, a logdet computation method, the order of polynomial to be used,
% and a number of samples. it then perfoms an experiment to try to maximize
% the log-likelihood function. 

%% Experiment 2
function [estimates, rho_vect] = log_likelyhood(fileName, rho_hidden, method, order, samp)
    %% evaluating input
    if nargin == 0
        fileName = 'kneser (17,8).txt';
        rho_hidden = -0.22;
        method = 'Rational';
        order = 14;
        samp = 20;
    elseif nargin == 5
        nk = split(fileName(9:end-5),',');
        fileNameValid = strcmp(fileName(1:8),'kneser (')&&(str2num(nk{1})==2*str2num(nk{2})+1)&&strcmp(fileName(end-4:end),').txt');
        if (order < 1 || samp < 1 ||...
                ~(strcmp(method,'Cholesky')||strcmp(method,'Chebyshev')||strcmp(method,'Rational')) ||...
                rho_hidden>1 || rho_hidden < -1 || ~fileNameValid)
            fprintf('invalid input\n');
            return
        else
            fprintf('using kneser graph and rho provided by user\n');
        end
    else
        return
    end
    
    %% Generating the Hidden Information Matrix Jrho0
    [n,k,A]=generate_adjacency_matrix(fileName);
    I = speye(size(A));
    mat_dim = nchoosek(n,k);
    degree = nchoosek(n-k,k);
    Jrho_hidden = A*rho_hidden + I*(degree+1); % making Jrho0 diagonally dominant by adding degree +1
   
    %% Random Samples of Jrho0
 	amd = symamd(Jrho_hidden); % The SYMAMD command uses the approximate minimum degree algorithm (a powerful graph-theoretic technique) to produce large blocks of zeros in the matrix.
    % a normally distibuted random sample of Jrho0 can be found using the
    % formula Lrho0'\randn(mat_dim,1). 
	fprintf('finding cholesky factorisation of Jrho0\n');
    Lrho0 = chol(Jrho_hidden(amd,amd)); 
    min_rho = -1;
    max_rho = 1;
    N = 100;
    step_size = (max_rho-min_rho)/N;
    rho_vect = min_rho:step_size:max_rho;
    estimates = zeros(size(rho_vect));
    fprintf('generating random samples\n');
    X = zeros(N,mat_dim);
    for rho_counter = 1:N % generate N random samples
        x = Lrho0'\randn(mat_dim,1);
        X(rho_counter,amd) = x'; % need to un-amd the x
        fprintf('found sample #%.f of %.f\n',rho_counter,N);
    end

    %% Sweeping over rho
    % Check values for rho in the range [-1,1]. the value that maximises
    % the log-likelihood function will be used as the estimate for rho0
    % maximum log-likelihood is found by maximising the equation
    % log_like = 0.5*N*logdetJrho - N*0.5*(Xi'*Jrho*Xi); 
    % and sweeping over rho, creating a new Jrho in each iteration
    for rho_counter=1:length(rho_vect)
        rho = rho_vect(rho_counter);
        Jrho = A*rho + I*(degree+1);
        
        % choose a method based on user input
        if strcmp(method,'Cholesky')
        	amd = symamd(Jrho); % The SYMAMD command uses the approximate minimum degree algorithm (a powerful graph-theoretic technique) to produce large blocks of zeros in the matrix.
            logdetJrho = full(2*sum(log(diag(chol(Jrho(amd,amd)))))); 
        elseif strcmp(method,'Rational')
            block_size = min(20,samp);
            logdetJrho = Rational(Jrho,order,samp,block_size);
        elseif strcmp(method,'Chebyshev')
            lmax = 5*max(max(Jrho));
            lmin = 0.85*min(2*diag(Jrho)'-sum(abs(Jrho)));
            logdetJrho = ChebLogDet(Jrho,samp,order,lmin,lmax);
        else
            return;
        end
                
%         log_like = 0.5*N*logdetJrho - N*0.5*(Xi'*Jrho*Xi);
        log_like=0.5*N*logdetJrho; 
        
        % -N*0.5*(Xi'*Jrho*Xi);
        for i=1:N
            Xi = X(i,:)';
            sample = -0.5*(Xi'*Jrho*Xi);
            log_like = log_like+sample;
        end
        estimates(rho_counter) = log_like;
        
        fprintf('checked rho = %.2f\n',rho);
    end
    
    
    %% Plotting the Results
    figure (1)
    plot(rho_vect, estimates)
    indexmax = find(max(estimates) == estimates);
    xmax = rho_vect(indexmax);
    ymax = estimates(indexmax);
    hold on
    plot(xmax,ymax,'ro')
    strmax = ['Maximum at \rho = ',num2str(xmax)];
    text(xmax,ymax,strmax,'HorizontalAlignment','right');
    hold off
    xlabel('\rho');
    ylabel('log-likelyhood estimate');
    title('Log-Likelyhood estimation for hidden \rho_0 = -0.22')
end