function a = constr_component_a(a, constrain_vector)
    [~, cols] = size(constrain_vector);
    
    for i = 1:cols
        if constrain_vector(1, i) == 1
            % For A
            temp = a(:, i);
            a(:, i) = 0;
            % a(i, i) = temp(i); % Ici il faut prendre en consideration c car si elle est xplicative de la variable lors tu passes a 0 
            temp = temp .* (i == (1:numel(temp)));
            a = [a, temp];
        end
    end
end

