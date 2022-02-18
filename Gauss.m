format compact;

% Calling the experiments with varoius generators and printing and showing
% both the outcomes in graph and table.
[a, b] = do_exper_gauss(9, @Generators.generate_1);
show_plot(a,b,"Eliminacja Gaussa dla danych 1.");

[c, d] = do_exper_gauss(9, @Generators.generate_2);
show_plot(c,d,"Eliminacja Gaussa dla danych 2.");

Generators.make_table(b,['Data 1'], a)
Generators.make_table(d,['Data 2'], c)


% Draws a plot for the given arguments
function show_plot(x, y, full_title)
    figure
    plot(x,y)
    
    grid on
    title(full_title, "FontSize", 18);
    xlabel('Rozmiar macierzy', "FontSize", 16);
    ylabel('Norma residum', "FontSize", 16);
end

% Gauss algorithm for the matrix
function out = gauss(matrix)
    %The starting point from where we try to find a bigger number in the
    %column
    len = length(matrix);
    for j = 1:len-2
        % Finding the biggest number in the column
        [number, position] = max(matrix(j, j:len -1));
        % Switching rows
        matrix([j position+j-1], :) = matrix([position+j-1 j],:);
        % Subtracting the current row from every row underneath
        for i = j+1:len-1
            m = matrix(i, j)/matrix(j, j);
            matrix(i, :) = matrix(i, :) - matrix(j, :)*m;
        end
    end
    out = matrix();
end

% Clears the matrix to a simple Ax = b formula
function out = Ax_b(matrix)
    N = length(matrix)-1;
    for i =N: -1 :1        
        for j = i -1 : -1 : 1  
            m = matrix(j,i)/matrix(i,i);
            matrix(j,:) = matrix(j, :) - matrix(i, :) * m;
        end
    end
    out = matrix;
end

% Count the values of the X's 
function out = count_x(matrix)
    N = length(matrix) - 1;
    out = zeros(N,1);
    out(N) = matrix(N, N+1)/matrix(N,N);
    for i = N -1: -1: 1
        out(i) = (matrix(i, N+1) - matrix(i,i+1:N)*out(i+1:N))/matrix(i,i);
    end
end


% Does the computaions for respective number of iterations, where with each
% iteration the matrix size grows twice.
function [size, residums] = do_exper_gauss(iter,  generator)
    start = 10;
    
    size = zeros(iter,1);
    residums = zeros(iter, 1);

    for i = 1:iter
        matrix = generator(start);
        calculated = gauss(matrix);
        x = count_x(calculated);

        % Calculating the residum norm
        sum = 0;
        for j = 1:start
            sum = sum + (calculated(j,1:start)*x(:) - calculated(j,start +1))^2;
        end
        residums(i) = sqrt(sum);

        size(i) = start;
        start = start*2;
    end
end


























