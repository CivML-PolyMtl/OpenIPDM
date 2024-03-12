function A = constr_component_a(a, constrain_vector, c)

c_bin = binary_matrix(c);
constrain_vector_bin = binary_matrix(constrain_vector);

temp = c_bin * a * constrain_vector_bin'; % Have to be test in other cases

[row_count, ~] = size(a);

temp_extended = zeros(row_count, size(temp, 2));
temp_extended(1:size(temp, 1), :) = temp;

disp(temp_extended);
A = [a, temp_extended];

row_indices = find(c == 1);
col_indices = find(constrain_vector == 1);

for i = 1:length(row_indices)
    for j = 1:length(col_indices)
        % Set the corresponding element of a to 0
        A(row_indices(i), col_indices(j)) = 0;
    end
end

end

