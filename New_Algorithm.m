function retVal = New_Algorithm(dataset, M, a, b)

[t,w] = lgwt(M,a,b);
dataset_size = size(dataset);

% logdet(B) ~ tr(log(B)) = 
%(1/N)*sum(i=1:N) (v'*[(B-I)*sum(j=1:m)(wj*(tj*B+(1-tj)*I)^-1]*v)
B = dataset;
I = speye(dataset_size);
v = ((rand(1,dataset_size(1))<.5)*2 - 1)'; % randmacher vector;

% pcg solves AX=b
% vi = (tj*B+(1-tj)*I)*b
N=7; % need to fully parameterise
% approx = (1/N) * (B-I).*sum(w./(t.*B + (1-t)*I));

% need to remove loops
for i=1:N
    grand_total = 0;
    total = 0;
    for j = 1:M
        wj = w(j);
        tj = t(j);
        hard_compute = (tj*B+(1-tj)*I)\v;
        calc = wj*hard_compute;
        total = total+calc;
    end
    BI_calc = (B-I)*calc;
    vBI_calc = v'*BI_calc;
    grand_total = grand_total + vBI_calc;
end
grand_total = grand_total/N;
retVal = grand_total;
end