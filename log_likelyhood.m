function log_likelyhood(fileName, rho_0)
    if nargin == 0
        fileName = 'kneser (13,6).txt';
        rho_0 = -0.22;
    elseif nargin == 2
        fprintf('using kneser graph and rho provided by user');
    else
        return
    end
    
    [n,k,A]=generate_adjacency_matrix(fileName);
    I = speye(size(A));
    N = nchoosek(n,k);
    degree = nchoosek(n-k,k);
    Jrho0 = A*rho_0 + I*(degree+1);
	amd = symamd(Jrho0); % The SYMAMD command uses the approximate minimum degree algorithm (a powerful graph-theoretic technique) to produce large blocks of zeros in the matrix.
	Lrho0 = chol(Jrho0(amd,amd));
%     Lrho0 = chol(Jrho0);
    
    
    min_rho = -1;
    max_rho = 1;
    no_steps = 20;
    step_size = (max_rho-min_rho)/no_steps;
    rho_vect = min_rho:step_size:max_rho;
    estimates = zeros(size(rho_vect));
    
    X = zeros(no_steps,N);
    for rho_counter = 1:no_steps
        x = Lrho0'\randn(N,1);
%         X(i,amd) = x(or amd here)'; % need to un-amd the x
        X(rho_counter,amd) = x'; % need to un-amd the x
    end
   
    for rho_counter=1:length(rho_vect)
        rho = rho_vect(rho_counter);
        Jrho = A*rho + I*(k+1);
        amd = symamd(Jrho); % The SYMAMD command uses the approximate minimum degree algorithm (a powerful graph-theoretic technique) to produce large blocks of zeros in the matrix.
        logdetJrho = full(2*sum(log(diag(chol(Jrho(amd,amd)))))); 
        log_like=0;
        for j=1:no_steps
            Xj = X(j,:)';
            sample = -0.5*(Xj'*Jrho*Xj)+0.5*logdetJrho;
            log_like = log_like+sample;
        end
        estimates(rho_counter) = log_like;
    end
    
    figure (1)
    plot(rho_vect, estimates)
    indexmax = find(max(estimates) == estimates);
    xmax = rho_vect(indexmax);
    ymax = estimates(indexmax);
    hold on
    plot(xmax,ymax,'ro')
    strmax = ['Maximum at rho = ',num2str(xmax)];
    text(xmax,ymax,strmax,'HorizontalAlignment','right');
    hold off
    xlabel('rho');
    ylabel('log-likelyhood');
    title('Log-Likelyhood estimation for hidden Rho')
end