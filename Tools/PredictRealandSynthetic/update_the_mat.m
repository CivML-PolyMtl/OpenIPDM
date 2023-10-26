function[flag] = update_the_mat(X)
    mat_type = 8;
    age_low= 1;
    age_high = 1000;
    lat_low = 50;
    lat_high = 52;
    djma_low = 1;%7000;
    djma_high = 1000000;%100000;
    trucks = 0;
    flag = false;
%                         && (X(1,7) < lat_high && X(1,7) > lat_low)...

%     long_left = -64.6;
%                     && X(1,8) < long_left...
%                         && X(1,5) == mat_type...
    if ~isempty(X)
        if ~isnan(X(1,13))
            if (X(1,6) < age_high && X(1,6) > age_low)...
                    && (X(1,13) < djma_high && X(1,13) > djma_low)...
                    && X(1,14) > trucks
                flag = true;
            end
        end
    end
end