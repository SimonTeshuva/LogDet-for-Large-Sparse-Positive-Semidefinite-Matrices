close all;
clear all;
clc;
% https://au.mathworks.com/help/matlab/ref/sprandsym.html
%need a far better method. this bricks my computer at somewhere between 1E4
%and 2E4

disp('starting');

dataset_size = 1e4;
epr = 10;
diag_dom_const = 1e-3;
tic
dataset = Sparse_Dataset(dataset_size, epr, diag_dom_const);
dataset_time = toc;
disp('dataset done')
n = 14;
Nm = 20;
M = 100;
a = 1e-3;
b = 1;
max_exact_size =3e4;
times = [];
vals = [];
percent_error = [];
fileID = fopen('output.txt', 'a+');

if dataset_size<3e4
    tic
    logdet_exact = full(2*sum(log(diag(chol(dataset)))));
    logdet_exact_time = toc;
end

for k = 1:1:floor(Nm/2)
    fprintf('checked %.f logdets\n', k);
    tic
    vals = [vals Rational(dataset, M, Nm, k)];
    times = [times toc];
    if dataset_size < max_exact_size
        percent_error = [percent_error abs(100*(1-vals(k)/logdet_exact))];
    end
end
tic
zero_blocks=Rational(dataset, M, Nm, Nm);
zero_blocks_time = toc;
zero_blocks_error = abs(100*(1-vals(k)/logdet_exact));

c = clock;
fprintf(fileID, 'results from %.f\t%.f\t%.f\t%.f\t%.f\t%.f on a dataset of size %.f\n', c, dataset_size);
fprintf(fileID,'k\t\ttimes\t\tlogdets\t\tpercentage error\t\tno_blocks\n');
if sum(percent_error) > 0
    fprintf(fileID,'exact\t%.5f\t%.2f\t\t0\t\t\t\t\t\tN/A\n',logdet_exact_time, logdet_exact);
    for counter = 1:length(times)
        fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',counter,times(counter),vals(counter), percent_error(counter), floor(Nm/counter));
    end
else
    for counter = 1:length(times)
        fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',counter,times(counter),vals(counter), 0, floor(Nm/counter));
    end
end
fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',Nm,zero_blocks_time, zero_blocks, zero_blocks_error, 1);
fprintf(fileID, '\n\n');

fclose(fileID);
% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory