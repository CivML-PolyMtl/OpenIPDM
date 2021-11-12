function [SumLL,grad_fun,InirilizedEx,InirilizedVar,X_ControlPoints,fx]=...
    AnalysisObjective(MdataEngy,InspectorsID,InspectorStructureIndex,...
    OptBoundsData,A,F, Q, x0, s2_X0,param,CurrentInspectorParam,...
    EngBiasData,CurrentInspectorID,Ncurve,DataTransformedSpace,...
    OptProcedure,ParallelComp,OptLevel,GlobalCondData,PriorKnowledgeChk,...
    MeanPrior, StdPrior,Nsigma,GradCompute,varargin)

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

SumLL=0;GPUCompute=0;kk=0;grad_fun=0;
if OptLevel==1
    EngBiasData(:,2)=param(2);
    CurrentInspectorID=[];
end

%% Compute the objective function value (sum of likelihoods)
if ~isempty(MdataEngy.YS)
    fx=zeros(size(MdataEngy.YS));
    if ParallelComp==2
        % Run Analyses on GPU
        GPUCompute=1;
            [fx,InirilizedEx,InirilizedVar,X_ControlPoints]=LogLikKF(InspectorsID,...
                EngBiasData,CurrentInspectorID,CurrentInspectorParam,...
                MdataEngy, A, F, Q(param),Ncurve,OptBoundsData,GlobalCondData,...
                OptProcedure,param,Nsigma,GPUCompute,kk,OptLevel,...
                StructuralAttributeAnalyses,Kernel_Length,Sigma_W0,...
                Kernel_Type,NumControlPoints,DoneStructuralAnalyses,...
                InitialEX,InitialVar,MdataEngy.NaNRemoved,TrueSpeedValues);
    elseif ParallelComp==1
        % Run Analyses in parallel on CPU 
        parfor kk=1:size(MdataEngy.YS,2)
            [fx(kk)]=LogLikKF(InspectorsID,EngBiasData,CurrentInspectorID,CurrentInspectorParam,...
                MdataEngy.AllInspectors,MdataEngy.YS(1,kk,:), A, F, Q(param), MdataEngy.RS(1,kk), MdataEngy.ReS(1,kk,:),...
                MdataEngy.InpecBiaseS(1,kk,:),MdataEngy.CurrentInspectorS(1,kk,:),MdataEngy.RUS(1,kk,:),...
                MdataEngy.InspBUS(1,kk),MdataEngy.ObsYearsS(1,kk,:),MdataEngy.yearlyS(1,kk,:),MdataEngy.InspectorLabelS(1,kk,:),...
                MdataEngy.StructureIndS(1,kk),MdataEngy.ElementIndS(1,kk),MdataEngy.init_x,Ncurve,OptBoundsData,...
                GlobalCondData,OptProcedure,param,Nsigma,GPUCompute,kk,OptLevel,StructuralAttributeAnalyses,Kernel_Length,Sigma_W0,Kernel_Type,NumControlPoints);
        end    
    else
        % Run Analyses on CPU 
        for kk=1:size(MdataEngy.YS,2)
            [fx(kk)]=LogLikKF(InspectorsID,EngBiasData,CurrentInspectorID,CurrentInspectorParam,...
                MdataEngy.AllInspectors,MdataEngy.YS(1,kk,:), A, F, Q(param), MdataEngy.RS(1,kk), MdataEngy.ReS(1,kk,:),...
                MdataEngy.InpecBiaseS(1,kk,:),MdataEngy.CurrentInspectorS(1,kk,:),MdataEngy.RUS(1,kk,:),...
                MdataEngy.InspBUS(1,kk),MdataEngy.ObsYearsS(1,kk,:),MdataEngy.yearlyS(1,kk,:),MdataEngy.InspectorLabelS(1,kk,:),...
                MdataEngy.StructureIndS(1,kk),MdataEngy.ElementIndS(1,kk),MdataEngy.init_x,Ncurve,OptBoundsData,...
                GlobalCondData,OptProcedure,param,Nsigma,GPUCompute,kk,OptLevel,StructuralAttributeAnalyses,Kernel_Length,Sigma_W0,Kernel_Type,NumControlPoints);
        end
    end
    % include the prior in the objective
if PriorKnowledgeChk==1
    stdp=StdPrior;
    mu=MeanPrior;
    s_ln=@(m,s) sqrt(log(1+(s/m)^2)); 
    m_ln=@(m,s_ln) log(m)-0.5*s_ln^2; 
    SL=s_ln(mu,stdp);
    ML=m_ln(mu,SL);
    if CurrentInspectorParam(1)~=0
        SumLL=sum(fx(:))+log(pdf('Lognormal',CurrentInspectorParam(1)-1,ML,SL));
    else
        SumLL=sum(fx(:))+log(pdf('Lognormal',param(2)-1,ML,SL));
    end
else
    SumLL=sum(fx(:));
end
end
%% compute the gradient of the objective function
if GradCompute==2
    SumLL=-SumLL;
elseif GradCompute==1
    SumLL=-SumLL;
    delta=[1E-6,1E-3,1E-3,1E-6,1E-6,1E-3,1E-3,1E-6];
    fx_TR=@(param) AnalysisObjective(MdataEngy,InspectorsID,...
        InspectorStructureIndex,OptBoundsData,A,F, Q, x0, s2_X0,param,...
        CurrentInspectorParam,EngBiasData,CurrentInspectorID,Ncurve,...
        DataTransformedSpace,OptProcedure,ParallelComp,OptLevel,...
        GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,Nsigma,0);
    Term1=@(x) fx_TR(x+max([[delta.*1E-4]' [delta.*abs(param)]'],[],2));
    Term2=@(x) fx_TR(x-max([[delta.*1E-4]' [delta.*abs(param)]'],[],2));
    Term3=@(x) max([[delta.*1E-4]' [delta.*abs(param)]'],[],2);
    grad_fct=@(x) (-Term1(x)+Term2(x))./(2*Term3(x));
    grad_fun=grad_fct(param);
end
end