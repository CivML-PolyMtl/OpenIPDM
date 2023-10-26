function[accept_time_series] = syn_fix(X)
    accept_time_series = true;
    step = false;
    for i=2:size(X,1)
        if step
            next_step = X(i,3)-X(i-1,3);
            if step ~= next_step || step ~=2
                accept_time_series=false;
                break
            end
            step = next_step;
        else
            step = X(i,3)-X(i-1,3);
        end
    end
end