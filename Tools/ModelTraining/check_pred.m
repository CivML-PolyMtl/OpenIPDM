function [nan_check, large_var_check, inf_variance] = check_pred(states)
%     [mz, Sz, ma, Sa] = tagi.extractStates(states);
    [mz, Sz] = tagi.extractStates(states);
    v2_indices = 2:2:size(mz{3,1});
    % nan_check = zeros(1,4);    
    nan_check = zeros(1,2);
    large_var_check = 0;
    inf_variance = 0;
    % nan check
    if any(isnan(mz{3,1}))
        nan_check(1,1) = 1;
    end
    % if any(isnan(ma{3,1}))
    %     nan_check(1,3) = 1;
    % end
    if any(isnan(Sz{3,1}))
        nan_check(1,2) = 1;
    end
    % if any(isnan(Sa{3,1}))
    %     nan_check(1,4) = 1;
    % end
    mv2a = mz{3,1}(v2_indices);
    Sv2a = Sz{3,1}(v2_indices);
    
    [mv2a, Sv2a, Cv2a] = act.expFun(mv2a, Sv2a, false);
    
    if Sv2a > 500
        large_var_check = 1;
    end
    if any(Sv2a == inf)
        inf_variance = 1;
    end
end

