% logdet exact formula is wrong
function Exp1(ex1_do)

if nargin==0
    ex1_do=[1,1,1,1];
elseif nargin==1
    fprintf('using user input\n');
else
    return
end
    

kappa=5;samp=20;order=14;
rs = [1e3,5e3,1e4,5e4,1e5,3e5,6e5,1e6,3e6,6e6,1e7];
rsra_chol=rs(rs<1e5);



if ex1_do(1)
    [crt,rrt,ckt,rkt,chol_rt,chol_kt]=runtime_comparison(rs,order,samp,kappa);
    k_vect=1:length(rs);
    n_vect=2*k_vect+1;
    for i=1:length(k_vect)
        ks(i)=nchoosek(n_vect(i),k_vect(i));
    end
    
    figure(1)
    plot(crt,log10(rs(1:length(crt))),rrt,log10(rs(1:length(rrt))),chol_rt,log10(rs(1:length(chol_rt))));
    ylabel('log10(Matrix Size)');
    xlabel('run time');
    title('Chebyshev,Rational,Cholesky on Random Matrix');
    legend('Chebyshev','Random','Cholesky');
    figure(2)
    plot(ckt,log10(ks(1:length(ckt))),rkt,log10(ks(1:length(rkt))),chol_kt,log10(ks(1:length(chol_kt))));
    ylabel('log10(Matrix Size)');
    xlabel('run time');
    title('Chebyshev,Rational,Cholesky on Kneser Based Matrix');
    legend('Chebyshev','Random','Cholesky');
end

if ex1_do(2)
    [rarr,rarc,rakr,rakc]=relative_accuracy_medium(rsra_chol,order,samp,kappa);
    k_vect=1:length(rs);
    n_vect=2*k_vect+1;
    for i=1:length(k_vect)
        ks(i)=nchoosek(n_vect(i),k_vect(i));
    end
    ks_ra_cheb=ks(1:length(rakr));
    
    figure(3)
    plot(log10(rsra_chol),rarr,log10(rsra_chol),rarc)
    xlabel('log10(Matrix Size)');
    ylabel('relative accuracy');
    title('Relative Accuracy on Random Matrix');
    legend('Chebyshev','Random');
    
    figure(4)
    plot(log10(ks_ra_cheb),rakr,log10(ks_ra_cheb),rakc)
    xlabel('log10(Matrix Size)');
    ylabel('relative accuracy');
    title('Relative Accuracy on Kneser Based Matrix');
    legend('Chebyshev','Random');
end

if ex1_do(3)
    k_vect=5:11;
    [rake_rat, rake_cheb]=relative_accuracy_large(k_vect,order,samp,kappa);
    n_vect=2*k_vect+1;
    for i=1:length(k_vect)
        ks(i)=nchoosek(n_vect(i),k_vect(i));
    end
    
    figure(5)
    plot(log10(ks),rake_rat,'bo-',log10(ks),rake_cheb,'r-');
    xlabel('log10(matrix size)');
    ylabel('relative accuracy');
    title('Relative Accuracy vs Exact Logdet from Kneser Based Matrix');
    legend('Rational', 'Chebysev');
end

if ex1_do(4)
    pertubations = [1e-4,1e-3,1e-2,1e-1,0.33,0.85,1,1.33,3,1e1,1e2,1e3,1e4];
    k=11;
    [rake_rat, rake_cheb]=bad_lambdas_vs_rational(k,order,samp,kappa,pertubations);    
    figure(6)
    plot(pertubations,rake_rat,'bo',pertubations,rake_cheb);
    xlabel('log10(pertubation)');
    ylabel('relative accuracy');
    title('Guessing lmin,lmax on Kneser Based Matrix');
    legend('Rational', 'Chebysev');
end
end

function [crt,rrt,ckt,rkt,chol_rt,chol_kt] = runtime_comparison(rs,order,samp,kappa)
diag_dom_const=1e-3;
ks=1:length(rs);
for i=1:max(ks)
    clear lmin lmax A Asize;
    Asize=rs(i);
    A=Sparse_Dataset(Asize,2*floor(log10(Asize)),diag_dom_const);
    if Asize<3e7
        lmin=min(2*diag(A)'-sum(abs(A)));lmax=max(max(A));
        tic;ChebLogDet(A,samp,order,lmin,lmax);crt(i)=toc;
    end
    if Asize<2e6
        block_size=min(20,order);
        parfor j=1
        end
        tic;Rational(A,order,samp,block_size);rrt(i)=toc;
    end
    if Asize<1e5
        tic;amd=symamd(A);full(2*sum(log(diag(chol(A(amd,amd))))));chol_rt(i)=toc;
    end
    
    clear lmin lmax A Asize;
    k=ks(i);n=2*k+1; Asize=nchoosek(n,k);
    fileName=kneser_file(k);
    A=generate_adjacency_matrix(fileName);
    lmin=-nchoosek(n-k-1,k); lmax=nchoosek(n-k,k);
    epsilon=(lmax-lmin)/(kappa+1);
    A=diag(A)+speye(size(A))*(epsilon-lmin);
    if Asize<3e7
        tic;ChebLogDet(A,samp,order,epsilon,lmax-lmin+epsilon);ckt(i)=toc;
    end
    if Asize<2e6
        parfor j=1
        end
        block_size=min(20,order);
        tic;Rational(A,order,samp,block_size);rkt(i)=toc;
    end
    if Asize<5e4
        tic;amd=symamd(A);full(2*sum(log(diag(chol(A(amd,amd))))));chol_kt(i)=toc;
    end
    fprintf('finished work on random matrix of size %.f and kneser graph with %.f nodes\n',rs(i),nchoosek(n,k))
end

end


function [rarr,rarc,rakr,rakc]=relative_accuracy_medium(rs,order,samp,kappa)
diag_dom_const=1e-3;
ks=1:length(rs);
for i=1:max(ks)
    
    clear lmin lmax A Asize;
    Asize=rs(i);
    A=Sparse_Dataset(Asize,2*floor(log10(Asize)),diag_dom_const);
    lmin=min(2*diag(A)'-sum(abs(A)));lmax=max(max(A));
    
    cld_r(i)=ChebLogDet(A,samp,order,lmin,lmax);
    block_size=min(20,order);
    parfor j=1
    end
    rld_r(i)=Rational(A,order,samp,block_size);
    
    amd=symamd(A);
    chol_ld_r(i)=full(2*sum(log(diag(chol(A(amd,amd))))));
    
    
    clear lmin lmax A Asize;
    k=ks(i);n=2*k+1; Asize=nchoosek(n,k);
    fileName=kneser_file(k);
    [n,k,A]=generate_adjacency_matrix(fileName);
    lmin=-nchoosek(n-k-1,k); lmax=nchoosek(n-k,k); % got a feeling this may not be legit
    epsilon=(lmax-lmin)/(kappa+1);
    A=diag(A)+speye(size(A))*(epsilon-lmin);
    
    cld_k(i)=ChebLogDet(A,samp,order,epsilon,lmax-lmin+epsilon);
    parfor j=1
    end
    block_size=min(20,order);
    rld_k(i)=Rational(A,order,samp,block_size);
    
    amd=symamd(A);
    chol_ld_k(i)=full(2*sum(log(diag(chol(A(amd,amd))))));
    
    fprintf('finished work on random matrix of size %.f and kneser graph with %.f nodes\n',rs(i),nchoosek(n,k))
end

rarr=1-(chol_ld_r-rld_r)./chol_ld_r;
rarc=1-(chol_ld_r-cld_r)./chol_ld_r;
rakr=1-(chol_ld_k-rld_k)./chol_ld_k;
rakc=1-(chol_ld_k-cld_k)./chol_ld_k;
end

function [rake_rat, rake_cheb]=relative_accuracy_large(k_vect,order,samp,kappa)

for i=1:length(k_vect)   
    clear lmin lmax A Asize;
    fileName=kneser_file(k_vect(i));
    [n,k,A]=generate_adjacency_matrix(fileName);
    lmin=-nchoosek(n-k-1,k); lmax=nchoosek(n-k,k); % got a feeling this may not be legit
    epsilon=(lmax-lmin)/(kappa+1);
    A=diag(A)+speye(size(A))*(epsilon-lmin);
    [ld,lmin2,lmax2]=ld_kneser(n,k,kappa);
    
    logdet_exact(i)=ld;
    
    cld_k(i)=ChebLogDet(A,samp,order,lmin2,lmax2);
    parfor j=1
    end
    block_size=min(20,order);
    rld_k(i)=Rational(A,order,samp,block_size);
    
    fprintf('finished work on kneser graph with %.f nodes\n',nchoosek(n,k))
    
end
    rake_rat = abs((logdet_exact-rld_k)./logdet_exact);
    rake_cheb = abs((logdet_exact-cld_k)./logdet_exact);
end

function [rake_rat, rake_cheb]=bad_lambdas_vs_rational(k,order,samp,kappa,pertubations)
clear lmin lmax A Asize;
n=2*k+1;
fileName=kneser_file(k);
A=generate_adjacency_matrix(fileName);
lmin=-nchoosek(n-k-1,k); lmax=nchoosek(n-k,k); % got a feeling this may not be legit
epsilon=(lmax-lmin)/(kappa+1);
A=diag(A)+speye(size(A))*(epsilon-lmin);
logdet_exact=ld_kneser(n,k,kappa);
for i=1:length(pertubations)
    lmin_pert = 1;
    lmax_pert = 1;
    if pertubations(i)<1
        lmin_pert = pertubations(i);
    elseif pertubations(i)>1
        lmax_pert = pertubations(i);
    end
        
    cld_k(i)=ChebLogDet(A,samp,order,lmin_pert*epsilon,lmax_pert*(lmax-lmin+epsilon));

    fprintf('finished work on kneser graph with %.f nodes\n',nchoosek(n,k));    
end

parfor j=1
end
block_size=min(20,order);
rld_k=Rational(A,order,samp,block_size);
rld_k=rld_k*ones(1,length(cld_k));

rake_rat = abs((logdet_exact-rld_k)./logdet_exact);
rake_cheb = abs((logdet_exact-cld_k)./logdet_exact);
end

function [logdet,lmin,lmax] = ld_kneser(n,k,kappa)
    evals = [];
    mults = [];
    for i = 0:k
        eval = ((-1)^i)*nchoosek(n-k-i, k-i); 
        if i == 0
            mult = 1;
        else 
            mult = nchoosek(n,i) - nchoosek(n, i-1);
        end
        evals = [evals eval];
        mults = [mults mult];
    end

	epsilon = (max(evals)-min(evals))/(kappa - 1);
    evals = evals - min(evals) + epsilon;
    
    logdet = sum(mults.*log(evals));
    lmin=min(evals);
    lmax=max(evals);
end