dt = 1.5;
A = [1, dt, (dt^2)/2; 0, 1, dt; 0, 0, 1];
constrain_vector = [0, 1, 0];

c = [1, 0, 0];
[row_count, col_count] = size(c);
c_ = [];

for i = 1:col_count
    if c(i) == 1
        % Créer une ligne de zéros
        row = zeros(1, col_count);
        % Mettre un 1 à la position correspondante
        row(i) = 1;
        % Ajouter cette ligne à la matrice de sortie
        c_ = [c_; row];
    end
end

c_*A*constrain_vector'



