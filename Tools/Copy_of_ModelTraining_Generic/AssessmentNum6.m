 OptimizationProceedure=1;

% initial parameter values
PARAM=[OptBoundsData(1,1) OptBoundsData(2,1) OptBoundsData(3,1) ...
    OptBoundsData(4,1) OptBoundsData(5,1) OptBoundsData(6,1)];

% Inspector uncertainty bounds
bounds={OptBoundsData(2,2:3)};

% Model parameter bounds
boundsQ={OptBoundsData(1,2:3) OptBoundsData(2,2:3) OptBoundsData(3,2:3)...
    OptBoundsData(4,2:3) OptBoundsData(5,2:3) OptBoundsData(6,2:3)};

if OperationIndex>=3
    % Inspector uncertainty bounds
    bounds={OptBoundsData(2,2:3) [-OptBoundsData(2,3) OptBoundsData(2,3)]};
end

OptimizationAlgorithmIndex=find(strcmp(app.SelectOptAlgorithm.Value,...
    app.SelectOptAlgorithm.Items));

% Initial optimization step
if isempty(ModelParameters)
    SigmaWSigmaVOneLevel();
    PARAM=Qparam;
    Stored_QParam = [];
else
    Qparam=cell2mat(ModelParameters);
    PARAM=Qparam;
    Stored_QParam = PARAM;
end

if get(app.LogParams,'value')==1 && isempty(ModelParameters)
    save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/2-PARAM%s%s.mat',pwd,Data_filename,date),'PARAM');
end
if get(app.LoglikelihoodPlot,'value')==1 && isempty(ModelParameters)
    cla(app.OutputPlot1);
    drawnow;
%     plot(app.OutputPlot1,datetime('now'),gather(fx_NR(end)),'o');
%     hold(app.OutputPlot1,'on');
%     grid(app.OutputPlot1,'on');
%     box(app.OutputPlot1,'on');
%     ylabel(app.OutputPlot1,'Logliklihood')
%     xlabel(app.OutputPlot1,'Time')
end
% Check the inspector parameters
if ~isempty(InspectorsParam)
    UpdatedInspectorsData=InspectorsParam{1};
    Stored_UpdatedInspectorData = UpdatedInspectorsData;
else 
    Stored_UpdatedInspectorData = [];
end
if app.DTransformedSpace.Value && ~isempty(app.TrueInspectorsData) && app.LogParams_2.Value
    cla(app.OutputPlot);
    xlim(app.OutputPlot,[1 7]);
    ylim(app.OutputPlot,[1 7]);
    plot(app.OutputPlot,[1,7],[1,7]);
    hold(app.OutputPlot,'on');
    plot(app.OutputPlot,UpdatedInspectorsData{1}(:,end),app.TrueInspectorsData{1}(:,1),'o');
    hold(app.OutputPlot,'off');
    grid(app.OutputPlot,'on')
    box(app.OutputPlot,'on');
    ylabel(app.OutputPlot,'Inspector True $\sigma(I_i)$','interpreter','latex')
    xlabel(app.OutputPlot,'Inspector Estimated $\sigma(I_i)$','interpreter','latex')
    drawnow;
elseif app.LogParams_2.Value
    histogram(app.OutputPlot,UpdatedInspectorsData{1}(:,end));
    xlabel(app.OutputPlot,'Histogram of Inspector Estimated Params $\sigma(I_i)$','interpreter','latex');
    drawnow;
end

if IncludeStructuralAtt
    OptLevel=3;
    boundKr=[];
    InitialEx=[];
    InitialVar=[];
    Kr_Param=KernelParameters{3}(:,1)';
    KRparam=zeros(2,1);
    StructuralAttributes=0;
    for vi=1:length(KernelParameters{3}(:,1))
        boundKr=[boundKr {KernelParameters{3}(vi,2:3)}];
    end
    if ~isempty(TrueStateData)
        TrueStateValuesLong=TrueStateData{2,1};
        TrueSpeedValues=cellfun(@(x) x(2,:), TrueStateValuesLong, 'UniformOutput', false);
        TrueSpeedValues(MdataEngy.RemovedData(:,1))=[];
        TrueSpeedValues=cellfun(@(x) x(1,1),TrueSpeedValues);
    else
        TrueSpeedValues=[];
    end
else
    RegressionModel=[];
    Stored_RegressionModel = [];
end

% Factoring Structural Attributes
% IncludeStructuralAtt
if ~isempty(RegressionModel) && IncludeStructuralAtt
    KRparam=zeros(1,1+length(RegressionModel.Kernel_l));
    InitialEx=RegressionModel.InirilizedEx;
    InitialVar=RegressionModel.InirilizedVar;
    KRparam(1)=RegressionModel.Sigma_W0;
    KernelParameters{1}=RegressionModel.KernelType;
    X_ControlPoints=RegressionModel.X_ControlPoints;
    KRparam(1,2:length(KRparam))=RegressionModel.Kernel_l;
    StructuralAttributes=MdataEngy.StrucAtt;
    Kr_Param=KRparam;
    Stored_RegressionModel = RegressionModel;
elseif ~IncludeStructuralAtt
    KRparam=zeros(2,1);
    KernelParameters{1}=0;
    KernelParameters{2}=0;
    InitialEx=[];
    InitialVar=[];
    Stored_RegressionModel = [];
end

    
% Full optimization
SigmaWSigmaVTwoLevels(); 



% Update the interface with the optimized parameter values
pause(0.1);
