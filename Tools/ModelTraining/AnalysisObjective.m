function [SumLL,grad_fun,InirilizedEx,InirilizedVar,X_ControlPoints,fx]=...
    AnalysisObjective(MdataEngy,InspectorsID,InspectorStructureIndex,...
    OptBoundsData,A,F, Q, x0, s2_X0,param,CurrentInspectorParam,...
    EngBiasData,CurrentInspectorID,Ncurve,OptProcedure,OptLevel,SpeedConstraints,...
    Nsigma,varargin)

args = varargin;
if ~isempty(args)
    StructuralAttributeAnalyses=args{1};
    Kernel_Length=args{2};
    Sigma_W0=args{3};
    Kernel_Type=args{4};
    NumControlPoints=args{5};
    TrueSpeedValues=args{8};
    if ~isempty(args{6})
        InitialEX=args{6};
        InitialVar=args{7};
        if sum(InitialEX~=0)>0
            % Evaluation
            DoneStructuralAnalyses=2;
        else
            % Post-training: computing: InitialEX, InitialVar
            DoneStructuralAnalyses=1;
        end
    else
        % Training
        DoneStructuralAnalyses=0;
        InitialEX=0;
        InitialVar=0;
    end
    if ~isfield('NaNRemoved',MdataEngy)
            MdataEngy.NaNRemoved=[];
    end
else
    StructuralAttributeAnalyses=0;
    Kernel_Length=0;
    Sigma_W0=0;
    Kernel_Type=0;
    NumControlPoints=0;
    DoneStructuralAnalyses=0;
    InitialEX=0;
    InitialVar=0;
    TrueSpeedValues=[];
    MdataEngy.NaNRemoved=[];
end

SumLL=0;kk=0;grad_fun=0;
if OptLevel==1
    EngBiasData(:,2)=param(2);
    CurrentInspectorID=[];
end

%% Compute the objective function value (sum of likelihoods)
if ~isempty(MdataEngy.YS)
    fx=zeros(size(MdataEngy.YS));
    % Run Analyses on GPU
    GPUCompute=1;
    [fx,InirilizedEx,InirilizedVar,X_ControlPoints]=LogLikKF(InspectorsID,...
        EngBiasData,CurrentInspectorID,CurrentInspectorParam,...
        MdataEngy, A, F, Q(param),Ncurve,OptBoundsData,SpeedConstraints,...
        OptProcedure,param,Nsigma,GPUCompute,kk,OptLevel,...
        StructuralAttributeAnalyses,Kernel_Length,Sigma_W0,...
        Kernel_Type,NumControlPoints,DoneStructuralAnalyses,...
        InitialEX,InitialVar,MdataEngy.NaNRemoved,TrueSpeedValues);
    SumLL=sum(fx(:));
end

end