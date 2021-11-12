InterventionType=1;
loglikelihood_Current=-10^6;
for IntervType=1:3
    if InterventionType~=0
        R_param=IntParam{1,IntervType}; % col: InterventionType
    else
        R_param=nan(1,6);
    end
    if isnan(R_param(1))
        R_param=NAParam{ImportanceInd,IntervType};
    end
    Q_r1=diag([R_param(1)^2 R_param(2)^2 R_param(3)^2]);                         % Local error (covariance)
    Q_r{1}=blkdiag(Q_r1,Q_r1);                                                       % Local error (full matrix)
    InterventionVar_Network{1}=VarIntParam{IntervType};                         % Prior effect (covariance)
    InterventionMu_Network{1}=ExIntParam{IntervType};
    % Default Int value
    if sum(isnan(InterventionVar_Network{1}(1,1)))>0
        InterventionVar_Network{1}=NAIntParam{ImportanceInd+2,IntervType};                                       % Prior effect (covariance)
        InterventionMu_Network{1}=NAIntParam{ImportanceInd,IntervType};
    end

    try
        [ExTry, VarKFTry, loglikelihood, ExsmoothTry, VsmoothTry, ~,...
                ~,~,~,~,~,LifeSpanIntTry]=KF_KS(y_Data,A,C,Q,Q_r,...
                Re,PriorParam,init_x,init_V,Ncurve,ConstrainedKF,InterventionCheck,...
                InterventionVector,InterventionMu_Network,InterventionVar_Network);
         if loglikelihood>loglikelihood_Current  
            Exsmooth=ExsmoothTry;
            Vsmooth=VsmoothTry;
            LifeSpanInt=LifeSpanIntTry;
            loglikelihood_Current=loglikelihood;
         end
    catch
        
    end
end