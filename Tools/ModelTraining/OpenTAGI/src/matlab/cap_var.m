function [Sz, Sa] = cap_var(Sz, Sa)
    v2_indices = 2:2:size(Sz{3,1});

    Sv2a = Sz{3,1}(v2_indices);
    Sv2a(Sv2a > 100) = 100;
    Sz{3,1}(v2_indices) = Sv2a;
    Sa{3,1}(v2_indices) = Sv2a;
end

