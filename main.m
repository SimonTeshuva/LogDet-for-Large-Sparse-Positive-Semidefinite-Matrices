close all;
clear all;
clc;
% https://au.mathworks.com/help/matlab/ref/sprandsym.html

disp('starting');

kappas = [1.01, 5, 1e5, 1e10];
orders = [2, 5, 14];
samples = [1, 5, 30];
lambda_low = [1, 0.85, 0.33, 0.01];
lambda_high= [1, 5, 25, 1e3]; 
averages = 1;

tic
for k = 11:11
    fileName = strcat('kneser (',num2str(2*k+1),',',num2str(k),').txt');
    Experiment3(fileName,kappas,orders,samples,lambda_low,lambda_high,averages);
end
experiment_time = toc;

% given lmin and lmax, the chebyshev method works well for some
% lmin_guess and lmax_guess when
% 0.47*lmin_guess, 0.47*lmax_guess < lmin_guess, lmax_guess < 7.5*lmin_guess, 7.5*lmax_guess
% ie; chebyshev method is accurate when condition number approximation is;
% .22*kappa < kappa_approx < 56*kappa
% outside those regions, rational method is more accurate

% diag_dom_const = 1e-3;
% dataset_size=1e4;
% epr=10;
% dataset=Sparse_Dataset(dataset_size,epr,diag_dom_const);
% n=14;
% m=100;
% logdet_cheb=Chebyshev(dataset,m,n,diag_dom_const);
% 
% 
% 
% menu()



% maxn = 12;
% for n = 10:maxn
%     tic; kneser = Generate_Kneser(2*n+1, n); time = toc;
%     fprintf('generating the (%.f, %.f) kneser graph took %.5f seconds\n', 2*n+1, n, time);
% end
% maxk = 9;
% logdets_kneser = zeros(1,maxk);
% times = zeros(1,maxk);
% for k = 1:maxk
%     n = 2*k+1;
%     fprintf('finding the reduced laplacian of the (%.f, %.f) kneser graph\n',n,k);
%     fileName = strcat('kneser (',num2str(n),',',num2str(k),').txt');
%     kneser = Laplacian(fileName);
%     fprintf('evaluating the logdet of the reduced laplacian of the (%.f, %.f) kneser graph''s adjacency matrix\n', n,k);
%     tic
%     logdets_kneser(k) = Rational(kneser, M, 36, 6);
%     times(k) = toc;
%     fprintf('the logdet of the reduced laplacian of the (%.f, %.f) kneser graph''s adjacency matrix\n is %.4f and took %.4f seconds to evaluate\n', n,k,logdets_kneser(k), times(k));
% end

% times = [];
% vals = [];
% percent_error = [];
% fileID = fopen('output.txt', 'a+');
% 
% if dataset_size<3e4
%     tic
%     logdet_exact = full(2*sum(log(diag(chol(dataset)))));
%     logdet_exact_time = toc;
% end
% 
% for k = 1:1:floor(Nm/2)
%     fprintf('checked %.f logdets\n', k);
%     tic
%     vals = [vals Rational(dataset, M, Nm, k)];
%     times = [times toc];
%     if dataset_size < max_exact_size
%         percent_error = [percent_error abs(100*(1-vals(k)/logdet_exact))];
%     end
% end
% tic
% zero_blocks=Rational(dataset, M, Nm, Nm);
% zero_blocks_time = toc;
% zero_blocks_error = abs(100*(1-vals(k)/logdet_exact));
% 
% c = clock;
% fprintf(fileID, 'results from %.f\t%.f\t%.f\t%.f\t%.f\t%.f on a dataset of size %.f\n', c, dataset_size);
% fprintf(fileID,'k\t\ttimes\t\tlogdets\t\tpercentage error\t\tno_blocks\n');
% if sum(percent_error) > 0
%     fprintf(fileID,'exact\t%.5f\t%.2f\t\t0\t\t\t\t\t\tN/A\n',logdet_exact_time, logdet_exact);
%     for counter = 1:length(times)
%         fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',counter,times(counter),vals(counter), percent_error(counter), floor(Nm/counter));
%     end
% else
%     for counter = 1:length(times)
%         fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',counter,times(counter),vals(counter), 0, floor(Nm/counter));
%     end
% end
% fprintf(fileID,'%.0f\t\t%.5f\t\t%.2f\t\t%.5f\t\t\t\t\t%.0f\n',Nm,zero_blocks_time, zero_blocks, zero_blocks_error, 1);
% fprintf(fileID, '\n\n');
% 
% fclose(fileID);
% rational approx https://www.mathworks.com/examples/matlab/community/22736-chebfun-guide-4-chebfun-and-approximation-theory