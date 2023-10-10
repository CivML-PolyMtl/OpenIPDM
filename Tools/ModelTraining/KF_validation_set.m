try
    [~, ~, ~, LogLikValidation,~,~] = kalman_filter(y_valid, A, C, Q, R_valid, Re_Valid...
        , init_x_valid, init_V_valid,InpecBiase_valid,CurrentInspectorObs_valid,...
        RU_valid,InspBU_valid,ObsYears_valid...
        ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
    break_ = false;
    if ~isempty(AnnModel_loaded) && RegressLoop == 1
        LogLikValidation_loaded = LogLikValidation;
    end
catch
    disp(' ')
    disp('Constraints Failiure (KF validation)')
    LogLikValidation=-99^10;
    loglik=-99^10;
    break_ = true;
end