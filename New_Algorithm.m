function retVal = New_Algorithm(dataset, M, a, b)

[t,w] = lgwt(M,a,b);

% logdet(B) ~ tr(log(B)) = 
%(1/N)*sum(i=1:N) (v'*[(B-I)*sum(j=1:m)(wj*(tj*B+(1-tj)*I)^-1]*v)
B = dataset;
I = speye(size(dataset));

v = ((rand(1,dataset_size(1))<.5)*2 - 1)'; % randmacher vector;

% pcg solves AX=b
% vi = (tj*B+(1-tj)*I)*b

approx = (B-I).*sum(w./(t.*x + (1-t)));

end