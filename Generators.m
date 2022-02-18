classdef Generators
    methods (Static)
        % Generates the matrix for the 1. exercise scenario.
        function out = generate_1(size)
                out = zeros(size, size+1);
                for i = 1:size 
                    for j = 1:size + 1
                        if i == j
                            out(i,j) = 8.0;
                        
                        elseif (i == j-1) || (i == j+1)
                            out(i,j) = 4.0;
                        
                        else
                            out(i,j) = 0.0;
                        end
                    end
                    out(i,size+1) = 4.0 + 0.3*i; 
                end
        end

        % Generates the matrix for the 2. exercise scenario.
        function out = generate_2(size)
            out = zeros(size, size+1);
            for i = 1:size 
                for j = 1:size + 1
                    out(i,j) = 6/(7*(i + j + 1));
                end
                if mod(i,2) == 0
                    out(i,size+1) = 1/(3*i);
                else
                    out(i,size+1) = 0;
                end
            end
        end
        % Generates a table of outcome
        function make_table(data, column, row)
            fig = uifigure;
            uit = uitable(fig);
            uit.Data = data;
            uit.ColumnName = column;
            uit.RowName = row;
        end
    end
end