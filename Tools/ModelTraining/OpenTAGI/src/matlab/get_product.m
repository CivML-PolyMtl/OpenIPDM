function result = get_product(z, sign_x, sign_y, ind)
    result = zeros(size(z));
    result(ind) = z(ind) .* sign_x(ind);
    result(~ind) = z(~ind) .* sign_y(~ind);
end
