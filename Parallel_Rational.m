function retVal = Parallel_Rational(dataset, M, N)
tmin = 0;
tmax = 1;
[t,w] = lgwt(M,tmin,tmax);
dataset_size = size(dataset);

% logdet(B) ~ tr(log(B)) =
%(1/N)*sum(i=1:N) (v'*[(B-I)*sum(j=1:m)(wj*(tj*B+(1-tj)*I)^-1]*v)
B = dataset;
I = speye(dataset_size);

% pcg solves AX=b
% vi = (tj*B+(1-tj)*I)*b
% approx = (1/N) * (B-I).*sum(w./(t.*B + (1-t)*I));

% instead of Sparse_Dataset, use a random matrix
% u = rand (number), z = u*u'
% feed z into new_algorithm. TESTING123

% \ method on sparse matrices may not be a good idea. investigate. try
% calling full() first

% need to remove loops
grand_total = 0;
% test sizes for N
parfor i=1:N
    v = ((rand(1,dataset_size(1))<.5)*2 - 1)'; % randmacher vector;
    total = 0;
    for j = 1:M
        wj = w(j);
        tj = t(j);
        % need to convert \ to pcg
        linear_system = (tj*B+(1-tj)*I)\v;
        total = total+wj*linear_system;
    end
    vBI_linsys = v'*(B-I)*total;
    grand_total = grand_total + vBI_linsys;
    fprintf('%.f\n', 100*i/N);
end
grand_total = grand_total/N;
retVal = grand_total;
end

% logdet appears to be inversly proportional to M and N