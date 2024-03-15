function [x, p] = relu_constrain_step(x, p, constrain_vector, lower_bounds, upper_bounds, smoother)
%{
    Apply ReLU constraint to the state vector and the covariance matrix.

    Args:
        x (np.ndarray): The state vector.
        p (np.ndarray): The covariance matrix of the state vector.
        constrain_vector (np.ndarray): The constraint vector.
        lower_bounds (np.ndarray): The lower bounds.
        upper_bounds (np.ndarray): The upper bounds.
        smoother(true or false): Do we have to perform the smoother 

    Returns:
        np.ndarray: The updated state vector.
        np.ndarray: The updated covariance matrix.
%}

disp('code')
for i = 1:size(constrain_vector, 2)
    H = zeros(1, size(p, 2));
    disp('Loop ok')
    if constrain_vector(i) == 1
        H(1, i) = 1;
        disp('if ok')
        [mu, variance, cov] = mixture_truncation(x, p, p(:, i), H, lower_bounds(i), upper_bounds(i));

        if smoother == false
            x = [x; mu];
            p = [p, cov; cov', variance];
        end

        if smoother == true
            x(i+size(constrain_vector, 2)-1) = mu;
            p(i+size(constrain_vector, 2)-1, :) = cov(:)';
            p(:, i+size(constrain_vector, 2)-1) = cov(:);
            p(i+size(constrain_vector, 2)-1, i+size(constrain_vector, 2)-1) = variance;
        end
    end
end

end

