function[accept_time_series] = count_obs(X)
    threshold = 10;
    if size(X,1) > threshold
        accept_time_series = true;
    else
        accept_time_series = false;
    end
end