function [X, x_mean_train, x_std_train, x_min, x_max] = TAGI_scale_data(X, normalize, min_max, min_, max_)
    if normalize
        [X, x_mean_train, x_std_train] = TAGI_scaler.fit_transform_norm(X);
        x_min = [];
        x_max = [];
    elseif min_max
        [X, x_min, x_max] = TAGI_scaler.fit_transform_min_max(X, min_, max_);
        x_mean_train = [];
        x_std_train = [];
    end
end