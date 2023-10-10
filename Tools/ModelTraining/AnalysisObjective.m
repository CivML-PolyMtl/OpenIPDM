function [SumLL,grad_fun,InirilizedEx,InirilizedVar,X_ControlPoints,fx, AnnModel_]=...
    AnalysisObjective(MdataEngy,InspectorsID,InspectorStructureIndex,...
    OptBoundsData,A,F, Q, x0, s2_X0,param,CurrentInspectorParam,...
    EngBiasData,CurrentInspectorID,Ncurve,...
    OptProcedure,OptLevel,GlobalCondData,...
    Nsigma,varargin)

args = varargin;
AnnModel = [];
ann_model_training = 0;

if ~isempty(args)
    StructuralAttributeAnalyses=args{1};
    Kernel_Length=args{2};
    Sigma_W0=args{3};
    Kernel_Type=args{4};
    NumControlPoints=args{5};
    TrueSpeedValues=args{8};
    % train/evaluate regression model
    if isempty(args{9}) % train/evaluate KR if the provided AnnModel is empty
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
    else % train/evaluate TAGI
        if isempty(args{9}.theta_nn)
            % train/evaluate
            AnnModel = args{9};
            if AnnModel.evaluate_mode
                ann_model_training = 0; % !!this might be wrong!!
            else
                ann_model_training = 1;
            end
        else
            AnnModel = args{9};
            if AnnModel.cont_training && ~AnnModel.evaluate_mode
                % cont training
                ann_model_training = 1;
            else
                % evaluate
                ann_model_training = 0;
            end
        end
        StructuralAttributeAnalyses=args{1};
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
else % no regression analyses
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
    AnnModel = [];
end

SumLL=0;GPUCompute=0;kk=0;grad_fun=0;
InirilizedEx=0;InirilizedVar=0;X_ControlPoints=0;
if OptLevel==1
    EngBiasData(:,end)=param(2);
    CurrentInspectorID=[];
end

%% Compute the objective function value (sum of likelihoods)
if ~isempty(MdataEngy.YS)
    % Run Analyses on GPU
    GPUCompute=1;
    [fx,InirilizedEx,InirilizedVar,X_ControlPoints, AnnModel_]=LogLikKF(InspectorsID,... % Question: need to assign the output of LogLikKF executed with TAGI here
        EngBiasData,CurrentInspectorID,CurrentInspectorParam,...
        MdataEngy, A, F, Q(param),Ncurve,OptBoundsData,GlobalCondData,...
        OptProcedure,param,Nsigma,GPUCompute,kk,OptLevel,...
        StructuralAttributeAnalyses,Kernel_Length,Sigma_W0,...
        Kernel_Type,NumControlPoints,DoneStructuralAnalyses,...
        InitialEX,InitialVar,MdataEngy.NaNRemoved,TrueSpeedValues,...
        ann_model_training, AnnModel);
        SumLL=sum(fx(:));
end

end