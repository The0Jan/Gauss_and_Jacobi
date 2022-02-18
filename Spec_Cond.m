% Performing the conditioning and counting the spectral radius for both
% kinds of matrixes of various sizes.

test = [10, 20, 40, 80, 160, 320, 640, 1280, 2560];
d1 = conditioning(test, @Generators.generate_1)
d2 = conditioning(test, @Generators.generate_2)

spec1 = spec_rad(test, @Generators.generate_1)
spec2 = spec_rad(test, @Generators.generate_2)

show_plot(d1,d2,test, "Miara wskaźnika uwarunkowania macierzy")

make_table([spec1; spec2], test, ['Dane 1' ; 'Dane 2'])


% Puts the output data into a tbale for further research.
function make_table(data, column, row)
    fig = uifigure;
    uit = uitable(fig);
    uit.Data = data;
    uit.ColumnName = column;
    uit.RowName = row;
end

% Performs Conditioning  given type of matrix
function out = conditioning(sizes, generator)
    conds = zeros(1,length(sizes));
    for i = 1:length(sizes)
        matrix = generator(sizes(i));
        matrix_y = matrix(:,1:end-1);
        inv_matrix = inv(matrix_y);
                
        conds(i) = norm(inv_matrix)*norm(matrix_y);
    end
    out = conds;
end

% Give the spectral radius
function out = spec_rad(sizes, generator)
    specs = zeros(1, length(sizes));

    for j = 1:length(sizes)
        matrix = generator(sizes(j));
        LU_matr = matrix(:,1:length(matrix)-1);
        B_matr = matrix(:,end);
 
        D_matr = zeros(length(LU_matr), length(LU_matr));
        % Preparing the D matrix, The L + U matrix and the B matrix
        for i = 1:length(LU_matr) 
            D_matr(i,i) = 1/LU_matr(i,i);
            B_matr(i) = B_matr(i)/LU_matr(i,i);
            LU_matr(i,i) = 0.0;
        end
        % Doing the -N(L+U) step
        M_matr = -1 * D_matr * LU_matr;
        specs(j) = max(abs(eig(M_matr)));
    end
    out = specs;
end


% Draws a plot for the given arguments
function show_plot(x1, x2, y, full_title)
    figure
    semilogy(y, x1)
    hold on;
    semilogy(y, x2)

    hold off
    legend('Macierz 1', 'Macierz 2', "Fontsize", 16);

    grid on
    title(full_title, "FontSize", 18);
    xlabel('Rozmiar macierzy', "FontSize", 16);
    ylabel('Wartość wskaźnika uwarunkowania', "FontSize", 16);
end