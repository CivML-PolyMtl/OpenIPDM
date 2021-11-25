function [loglik,InirilizedEx,InirilizedVar,X_ControlPoints]=LogLikKF(InspectorsID,...
    EngBiasData,CurrentInspectorID,CurrentInspectorParam,MdataEngy,...
    A, C, Q,Ncurve,OptBoundsData,GlobalCondData,OptProcedure,param,...
    Nsigma,GPUCompute,kk,OptLevel,varargin)

InspectorsObs=MdataEngy.AllInspectors;
y=MdataEngy.YS;
R=MdataEngy.RS; 
Re=MdataEngy.ReS;
InpecBiase=MdataEngy.InpecBiaseS;
OptmInsp=MdataEngy.CurrentInspectorS;
RU=MdataEngy.RUS;
InspBU=MdataEngy.InspBUS;
ObsYears=MdataEngy.ObsYearsS;
yearly=MdataEngy.yearlyS;
InspectorLabel=MdataEngy.InspectorLabelS;
StructureInd=MdataEngy.StructureIndS;
ElementInd=MdataEngy.ElementIndS;
IniObsAvg=MdataEngy.init_x;

args = varargin;
% Initial values
StructuralAttributeAnalyses=0;
InirilizedEx=0;
InirilizedVar=0;
X_ControlPoints=0;
TrueSpeedValues=[];
if ~isempty(args)
    if sum(sum(sum(args{1}(1,:)~=0)))>0
        StrucAtt=args{1};
        Kernel_l=args{2};
        Sigma_W0=args{3};
        KernelType=args{4};
        Ncp=args{5};
        NanValues=args{9};
        TrueSpeedValues=args{10};
        MultiPass=10;
        % args{6}: Training=0, Post-training=1, Evaluation=2;
        if args{6}==0
            StructuralAttributeAnalyses=1;
        elseif args{6}==2
            % Evaluation
            MultiPass=2;
            InirilizedEx=args{7};
            InirilizedVar=args{8};
        end
    end
end

if GPUCompute
    for i=1:length(InspectorsID)
        ReReshape=reshape(Re,[size(Re,2),size(Re,3)]);
        ReReshape(find(InspectorsObs{i}))=EngBiasData(i,end).^2;
        Re(1,:,:)=ReReshape;
        
        % update inspector bias
        InpecBiaseReshape=reshape(InpecBiase,[size(InpecBiase,2),size(InpecBiase,3)]);
        InpecBiaseReshape(find(InspectorsObs{i}))=EngBiasData(i,2);
        InpecBiase(1,:,:)=InpecBiaseReshape;
    end
    CurrentInspectorObs=zeros(size(InspectorsObs{1}),'gpuArray');
    if ~isempty(CurrentInspectorID)
        RUReshape=reshape(RU,[size(RU,2),size(RU,3)]);
        CurrentInspectorIndex=find(CurrentInspectorID==InspectorsID);
        CurrentInspectorObs=InspectorsObs{CurrentInspectorIndex};
        RUReshape(find(CurrentInspectorObs))=(CurrentInspectorParam(1)).^2;
        RU(1,:,:)=RUReshape;
        
        if length(CurrentInspectorParam) ~= 1
            % update inspector bias
            InspBUReshape=reshape(InspBU,[size(InspBU,2),size(InspBU,3)]);
            InspBUReshape(find(CurrentInspectorObs))=(CurrentInspectorParam(2));
            InspBU(1,:,:)=InspBUReshape;
        end

    end
    
    Q=gpuArray(Q);
    A=gpuArray(A);
    C=gpuArray(C);
else
    for i=1:length(InspectorsID)
        ReReshape=reshape(Re,[size(Re,2),size(Re,3)]);
        ReReshape(find(InspectorsObs{i}(kk,:)))=EngBiasData(i,end).^2;
        Re(1,:,:)=ReReshape;
    end
    CurrentInspectorObs=zeros(size(InspectorsObs{1}(kk,:)));
    if ~isempty(CurrentInspectorID)
        RUReshape=reshape(RU,[size(RU,2),size(RU,3)]);
        CurrentInspectorIndex=find(CurrentInspectorID==InspectorsID);
        CurrentInspectorObs=InspectorsObs{CurrentInspectorIndex}(kk,:);
        RUReshape(find(CurrentInspectorObs))=(CurrentInspectorParam(1)).^2;
        RU(1,:,:)=RUReshape;
    end
end
if args{6}>=1
    % Evaluation
    KernelRegression_StructuralAtt1;
else
    % Include Structural Attributes
    if StructuralAttributeAnalyses
        % Trainingg
        KernelRegression_StructuralAtt1;
    else
        KalmanFilterSmoother_Parameters;
    end
end
end
