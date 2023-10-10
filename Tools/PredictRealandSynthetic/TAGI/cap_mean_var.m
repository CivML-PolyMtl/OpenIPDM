function [mz, Sz, ma, Sa] = cap_mean_var(mz, Sz, ma, Sa)
    v2_indices = 2:2:size(Sz{3,1});
    
    mv2a = mz{3,1}(v2_indices);
    Sv2a = Sz{3,1}(v2_indices);
    Sv2a(Sv2a > 10) = 10;
    mv2a(mv2a > 5) = 5;

    Sz{3,1}(v2_indices) = Sv2a;
    Sa{3,1}(v2_indices) = Sv2a;
    mz{3,1}(v2_indices) = mv2a;
    ma{3,1}(v2_indices) = mv2a;
end

