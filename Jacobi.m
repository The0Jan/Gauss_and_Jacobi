format compact;
% Perform the Jacobi method for various matrix sizes and matrix kinds.

iterations = [10 100 1000 5000 10000 50000];

%do_exper_jacobi(4,iterations,@Generators.generate_1,'Jacobi Algorithm for data 1.');
[matrix_size_1, residium_norm_1] = do_exper_jacobi(8,iterations,@Generators.generate_1,'Metoda Jacobiego dla danych 1.');
%do_exper_jacobi(4,iterations,@Generators.generate_1,'Jacobi Algorithm for data 2.');
[matrix_size_2, residium_norm_2] = do_exper_jacobi(8,iterations ,@Generators.generate_2,'Metoda Jacobiego dla danych 2.');

% Make a table out of the outcome
%Generators.make_table(residium_norm_1,iterations,matrix_size_1)
%Generators.make_table(residium_norm_2,iterations,matrix_size_2)





% Shows the plot done experiment for the Jacobi method
function show_plot_j(x, y, j_iter, full_title)
    figure
    legends = strings(length(j_iter),1);
    for i= 1:length(j_iter)
        plot(x,y(:,i));
        legends(i) = num2str(j_iter(i));
        hold on
    end
    grid on
    lgd = legend( legends(:), "FontSize", 14);
    
    title(full_title, "FontSize", 18);
    xlabel('Rozmiar macierzy', "FontSize", 16);
    ylabel('Norma residum', "FontSize", 16);
end


% Do a series of tests using the Jacobi method.
% Takes as input: N for the maximal matrix size(10 * 2^N), an array of
% numbers of iterations on the Jacobi method, and a matrix generator.
function [sizes residums] = do_exper_jacobi(iter, j_iter, generator, title)
    start = 10;

    sizes = zeros(iter,1);
    residums = zeros(iter, length(j_iter));

    for i = 1:iter
        for j = 1:length(j_iter)
            matrix = generator(start);
            residums(i, j) = jacobi( matrix, j_iter(j));

        end
        sizes(i) = start;
        start = start*2;
    end
    show_plot_j(sizes, residums, j_iter, title);

end


% A function that performs the Jacobi method on a given matrix, 
% with a given number of iterations.
function out = jacobi(matrix, iterations)
    alfa_1 = 2^-30;
    alfa_2 = 2^-30;
    
    LU_matr = matrix(:,1:length(matrix)-1);
    B_matr = matrix(:,end);
 
    D_matr = zeros(length(LU_matr), length(LU_matr));
    % Preparing the D matrix, The L + U matrix and the B matrix
    for i = 1:length(LU_matr) 
        D_matr(i,i) = 1/LU_matr(i,i);
        B_matr(i) = B_matr(i)/LU_matr(i,i);
        LU_matr(i,i) = 0.0;
    end
    % Doing the M =-N(L+U) step
    M_matr = -1 * D_matr * LU_matr;
    
    % Iterating over the X's
    X_mat_pre = zeros(length(LU_matr), 1);
    X_mat_after = zeros(length(LU_matr), 1);
    for round = 1:iterations
        for i = 1:length(LU_matr)
            X_mat_after(i) = B_matr(i) + M_matr(i,:)*X_mat_pre(:);
        end

        % Stop tests
        X_delta = X_mat_after-X_mat_pre;
        % The X delta stop test
        if norm(X_delta) < alfa_1
            % Counting the norm 
            norm_r = 0;
            for i = 1:length(matrix) -1
                norm_r = norm_r + (matrix(i,1:end-1)*X_mat_after(:) - matrix(i, end))^2;
            end
            norm_r = sqrt(norm_r);

            % The norm value stop test
            if norm_r < alfa_2
                out = norm_r;
                return ;
            else
                alfa_1 = alfa_1/2;
            end
        end
        X_mat_pre = X_mat_after;
        % The 'NaN' value stop test
        if sum(isnan(X_mat_pre)) 
            out = NaN;
            return
        end
    end
    out = 0;
    % Counting the residuum norm
    for i = 1:length(matrix) -1
        out = out + (matrix(i,1:end-1)*X_mat_pre(:) - matrix(i, end))^2;
    end
    out = sqrt(out);
end

