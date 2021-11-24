function [loglik,MAcc,AccVal,AccSTD,PWithinCI,PWithinCITrue,Erbar,Erbarbar,...
    BoundViolation,XtrueValue,XValue,XbtrueTr,XbbtrueTr,x,s_Xsmooth,V,Erbias,...
    Erbarbias, AccValBias]=KFsKF(app,y,...
    A, C, Q, R, Re, init_x, init_V,InpecBiase,InpecBiaseTrue,OptmInsp,RU,InspBU,ObsYears,...
    yearly,InspectorLabel,StructureInd,ElementInd,Pn,ReTrue,QTrue,...
    SynDatabaseState,RegressionModel,AllAtt,TableOfParameters,varargin)

%% Initilize Framework
InitilizeModel();

%% Run Kalman Filter
[x, V, ~, loglik,~,~,~,~,~,~] = kalman_filter(y, A, C, Q, R, Re, init_x, init_V,...
    InpecBiase,OptmInsp,RU,InspBU,ObsYears,Pn, OptBoundsData,...
    GlobalCondData(3,1),GlobalCondData,Nsigma);

%% Run Kalman Smoother
[Exsmooth,Vsmooth,s_Xsmooth,Status]=kalman_smoother(x,V,A,Q,Tdb,RunSmoother);

if isempty(RegressionModel)
ConditionVal=Exsmooth(1,2);
if isnan(Exsmooth(1,2))
    Exsmooth(1,2)
end
%% Refining initial state estimate (condition)
MaxCondition=100;
[Mtrv]=RevSpaceTransform(Pn,ConditionVal);
DifferenceObs=MaxCondition-Mtrv;
    if isempty(app)
        init_x(2)=TableOfParameters{2,5}*DifferenceObs;
        init_x(3)=TableOfParameters{2,6}*init_x(2);
        init_V(1,1,:)=max(TableOfParameters{2,2}^2,Re(2));
        init_V(2,2,:)=(TableOfParameters{2,3})^2.*(DifferenceObs)+...
            (TableOfParameters{2,7}).^2;
        init_V(3,3,:)=(TableOfParameters{2,4}.^2);
    else
        init_x(2)=TableOfParameters(1,5)*DifferenceObs;
        init_x(3)=TableOfParameters(1,6)*init_x(2);
        init_V(1,1,:)=max(TableOfParameters(1,2).^2,Re(2));
        init_V(2,2,:)=(TableOfParameters(1,3)).^2.*(DifferenceObs)+...
            (TableOfParameters(1,7)).^2;
        init_V(3,3,:)=(TableOfParameters(1,4).^2);
    end
    %% Run Kalman Filter
    [x, V, ~, loglik,~,~,~,~,~,~] = kalman_filter(y, A, C, Q,...
        R, Re, init_x, init_V,InpecBiase,OptmInsp,RU,InspBU,ObsYears,Pn,...
        OptBoundsData,GlobalCondData(3,1),GlobalCondData,Nsigma);

    %% Run Kalman Smoother
    [Exsmooth,Vsmooth,s_Xsmooth,Status]=kalman_smoother(x,V,A,Q,Tdb,RunSmoother);
end


MVv=s_Xsmooth;
x(:,1:Tdb)=Exsmooth(:,1:Tdb);

%% Run Kalman Filter - Perform Analyses with true parameters
if ~isempty(ReTrue)
    % Run Kalman Filter
    [xTrue, VTrue, ~, ~,~,~,~,~,~,~] = kalman_filter(y, A, C, QTrue, R,...
        ReTrue, init_xTrue, init_VTrue,InpecBiaseTrue,OptmInsp,...
        RU,InspBU,ObsYears,Pn,OptBoundsData,GlobalCondData(3,1),...
        GlobalCondData,Nsigma);
    % smoother pass
    [ExsmoothTrue,VsmoothTrue,s_XsmoothTrue,Status]=kalman_smoother(xTrue,...
        VTrue,A,QTrue,Tdb,RunSmoother);
    MVvTrue=s_XsmoothTrue;
end

%% Orgnize and Plot Output
rng(6)
OrgnizeOutputData();
PlotOutputData();

end