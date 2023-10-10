classdef TAGI_scaler
    methods(Static)
        function [x_mean_train, x_std_train] = fit_norm(X)
            x_mean_train = mean(X, 'omitnan');
            x_std_train = std(X, 'omitnan');
        end
        function x_scaled = transform_norm(x, x_mean, x_std)
            x_scaled = (x - x_mean)./x_std;
        end
        function [x_scaled, x_mean_train, x_std_train] = fit_transform_norm(X)
            x_mean_train = mean(X,1,'omitnan');
            x_std_train = std(X,1,'omitnan');
            x_scaled = (X - x_mean_train)./x_std_train;
        end
        function [x, x_var] = inverse_transform_norm(x_scaled, x_std_scaled, x_train_mean, x_train_std)
            x = x_scaled .* x_train_std + x_train_mean;
            if ~isempty(x_std_scaled)
                x_var = (x_std_scaled.^2) .* (x_train_std.^2);
            else
                x_var  = [];
            end
        end
        function x_scaled = transform_min_max(x, x_min, x_max, min_, max_)
            x_std = (x - x_min) ./ (x_max - x_min);
            x_scaled = x_std .* (max_ - min_) + min_;
        end
        function [x_min, x_max] = fit_min_max(x)
            x_min = min(x);
            x_max = max(x);
        end
        function [x_scaled, x_min, x_max] = fit_transform_min_max(x, min_, max_)
            x_min = min(x,[],1);
            x_max = max(x,[],1);
            x_std = (x - x_min) ./ (x_max - x_min);
            x_scaled = x_std .* (max_ - min_) + min_;
        end
        function x = inverse_transform_min_max(x_scaled, x_min, x_max, min_, max_)
            x_std = (x_scaled - min_) ./ (max_ - min_);
            x = x_std .* (x_max - x_min) + x_min;
        end
        function [x_scaled, x_min, x_max, var_factor] = fit_transform_min_max_with_var(x, min_, max_)
            x_min = min(x);
            x_max = max(x);
            x_std = (x - min(x)) ./ (max(x) - min(x));
            x_scaled = x_std .* (max_ - min_) + min_;
            x_var_train = var(x, 'omitnan');
            x_var_scaled = var(x_scaled, 'omitnan');
            var_factor = x_var_scaled/x_var_train;
        end
        function [x, x_var] = inverse_transform_min_max_with_var(x_scaled, x_var_scaled, x_min, x_max, min_, max_, var_factor)
            x_std = (x_scaled - min_) ./ (max_ - min_);
            x = x_std .* (x_max - x_min) + x_min;
            x_var = x_var_scaled/var_factor;
        end
    end
end
