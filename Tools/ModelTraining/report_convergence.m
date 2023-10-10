if Converge
%     if stalled
%         disp('Stall limit reached')
%     else
%         fprintf('Regress loop converged...\n');
%     end
    fprintf('Regress loop converged...\n');
    fprintf('Best performance on regress loop %d\n', Best_RegressLoop);
    end_time = toc(time_start);
    if total_time_TAGI == 0
        total_time_TAGI = end_time;
    else
        total_time_TAGI = total_time_TAGI + end_time;
    end
    RegressLoop=MultiPass;
else
    fprintf('Max multipass iterations exceeded\n');
    fprintf('Best performance on regress loop %d\n', Best_RegressLoop);
    end_time = toc(time_start);
    if total_time_TAGI == 0
        total_time_TAGI = end_time;
    else
        total_time_TAGI = total_time_TAGI + end_time;
    end
    Converge=1;
end