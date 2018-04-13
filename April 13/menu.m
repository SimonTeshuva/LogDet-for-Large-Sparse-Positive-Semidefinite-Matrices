%% Menu for FYP
% Author: Simon Tesuva
% Date Last Modified: 5/4/2018

% this function operates a menu system to allow the user to run various
% functions that were developed in this project. functions that can be run
% include: Sparse_Dataset(...), Generate_Kneser(...), Laplacian(...),
% Chebyshev(...) and Rational(...)

%% Code
function menu()
diag_dom_const = 1e-3;
dataset = [];
reduced_laplacian = [];

choice = -1;
while choice~=6
    fprintf('1 = Generate Diagonally Dominant Large Sparse Positive Semidefinite Matrix\n');
    fprintf('2 = Generate Kneser Graph (2k+1, k)\n');
    fprintf('3 = Generate Reduced Laplacian of a (2k+1, k) Kneser Graph\n');
    fprintf('4 = Use the Chebyshev Method\n');
    fprintf('5 = Use the Rational Function Method\n');
    fprintf('6 = quit\n');
    choice = (input('what to do?\n'));
    switch(choice)
        case 1
            dataset_size = (input('how large should the matrix be?\n'));
            epr = (input('how many elements per row should the matrix have'));
            tic
            dataset = Sparse_Dataset(dataset_size, epr, diag_dom_const);
            dataset_time = toc;
            disp('dataset done');
        case 2
            k = (input('what should k be?\n'));
            kneser = Generate_Kneser(2*k+1, k);
        case 3
            maxk = 9;
            toprint = strcat('which kneser?\nchoose a value for k from 1 to  ',num2str(maxk),'\n');
            choice2 = (input(toprint));
            if choice2 >= 1 && choice2 <= maxk
                k = choice2;
                kneser_name = strcat('kneser (',num2str(2*k+1),',',num2str(k),').txt');
                reduced_laplacian = Laplacian(kneser_name);
            else
                fprintf('this kneser graph has not yet been made.\n');
            end
        case 4
            choice2=(input('1=diagonaly dominant matrix\n2=reduced laplacian'));
            matrix = [];
            if choice2==1
                if ~isempty(dataset)
                    matrix=dataset;
                else
                    fprintf('A matrix has not yet been made.\n');
                end
            elseif choice2==2
                if ~isempty(reduced_laplacian)
                    matrix=reduced_laplacian;
                else
                    fprintf('A matrix has not yet been made.\n');
                end
            else
                fprintf('invalid choice\n');
            end
            if ~isempty(matrix)
                m = (round(input('pick a value for m, the sampling number\n')));
                n = (round(input('pick a value for n, the degree of the chebyshev polynomial\n')));
                tic
                logdet_chebyshev = Chebyshev(matrix,m,n,diag_dom_const);
                time=toc;
                if choice2==1
                    matrix_type='diagonally dominant';
                    matrix_size = dataset_size;
                elseif choice2==2
                    matrix_type='reduced laplacian';
                    matrix_size = nchoosek(n,k);
                else
                    matrix_type='invalid';
                end
                if ~strcmp(matrix_type,'invalid')
                    fprintf('the logdet of the %s matrix of size %.f by %.f is \n%.4f and took %.4f seconds to calculate\n',matrix_type,matrix_size,matrix_size,logdet_chebyshev,time);
                end
            end
        case 5
            choice2=(input('1=diagonaly dominant matrix\n2=reduced laplacian\n'));
            matrix = [];
            if choice2==1
                if ~isempty(dataset)
                    matrix=dataset;
                else
                    fprintf('A matrix has not yet been made.\n');
                end
            elseif choice2==2
                if ~isempty(reduced_laplacian)
                    matrix=reduced_laplacian;
                else
                    fprintf('A matrix has not yet been made.\n');
                end
            else
                fprintf('invalid choice\n');
            end
            if ~isempty(matrix)
                N = (round(input('pick a value for N, the sampling number\n')));
                M = (round(input('pick a value for M, the order of the legandre gradiant rational function\n')));
                k = (round(input('pick a value for k, the block size\n')));
                tic
                logdet_rational = Rational(matrix,M,N,k);
                time = toc;
                if choice2==1
                    matrix_type='diagonally dominant';
                    matrix_size = dataset_size;
                elseif choice2==2
                    matrix_type='reduced laplacian';
                    matrix_size = nchoosek(2*k+1,k)-1;
                else
                    matrix_type='invalid';
                end
                if ~strcmp(matrix_type,'invalid')
                    fprintf('the logdet of the %s matrix of size %.f by %.f is \n%.4f and took %.4f seconds to calculate\n',matrix_type,matrix_size,matrix_size,logdet_rational,time);
                end
            end
        case 6
            break
        otherwise
            fprintf('invalid choice\n');
    end
end

end