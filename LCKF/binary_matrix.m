function [output_matrix] = binary_matrix(matrix)

[~, col_count] = size(matrix);
output_matrix = [];

for i = 1:col_count
    if matrix(i) == 1
        % Créer une ligne de zéros
        row = zeros(1, col_count);
        % Mettre un 1 à la position correspondante
        row(i) = 1;
        % Ajouter cette ligne à la matrice de sortie
        output_matrix = [output_matrix; row];
    end
end
end

