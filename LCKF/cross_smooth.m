function [cross_cov] = cross_smooth(x_pred,p_pred,constrain_vector, cross_cov, lower_bounds, upper_bounds)

for i = 1:size(constrain_vector, 2)
    H = zeros(1, size(p_pred, 2));
    if constrain_vector(i) == 1
        H(1, i) = 1;

        [~, ~, cov] = mixture_truncation(x_pred, p_pred, cross_cov(:, i), H, lower_bounds(i), upper_bounds(i));
        cross_cov = [cross_cov, cov];
        disp(cross_cov)
    end
end
   
end

