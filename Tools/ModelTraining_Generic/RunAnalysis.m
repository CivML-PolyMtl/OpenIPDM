
global OptBoundsData GlobalCondData TrueInspectorsData KernelParameters MdataEngy AnalysesStartYear TestPercent
global UpdatedTrueStateData TrueStateAnalysisType TSStartingYear TSAnalysisDuration PriorKnowledgeChk MeanPrior StdPrior TrueStateData

OptimizationProceedure=find(strcmp(app.SelectOptimization.Value,app.SelectOptimization.Items));
OperationIndex=find(strcmp(app.SelectOptimization.Value,app.SelectOptimization.Items));
% Analysis Compute options
ComputeOptions=app.SelectOptAlgorithm.Value;
if find(strcmp(app.SelectOptAlgorithm.Value,app.SelectOptAlgorithm.Items))==2
    ParComp=1;
else
    ParComp=2;
end
GradCompute=0;  % 1 for trust region optimization
InspStrucIndex=[];%FindInspStrucIndex(VIData{1},UpdatedInspectorsData{1}(:,1));

if AnalysisSpace==1
    NTr=[];
end
%% Kalman Filter
dt=1;
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q=@(param) param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;(dt^4)/8 (dt^3)/3 (dt^2)/2;(dt^3)/6 (dt^2)/2 dt];
F=[1,0,0]; 
x0=zeros(3,1);
s2_X0=zeros(3);
UpdatedInspectorsData{1}(:,2)=OptBoundsData(2,1);
if OperationIndex>1
    IncludeStructuralAtt=1;
    NumStructuralAtt=KernelParameters{4};
else
    IncludeStructuralAtt=0;
    NumStructuralAtt=[];
end
if isempty(MdataEngy)
    ThisYear=date;
    ThisYear=str2double(ThisYear(end-3:end));
    TrainingWindow=ThisYear-AnalysesStartYear;
    [MdataEngy]=OrgnizeData(app.VIData{1},90E10,[0 0],UpdatedInspectorsData{1},...
        DataTransformedSpace,NTr,OptBoundsData,TrainingWindow,ParComp,...
        IncludeStructuralAtt,NumStructuralAtt,TestPercent);
end
if IncludeStructuralAtt && isempty(RegressionModel)
    if KernelParameters{5}
        RangeValues=gather((range(MdataEngy.StrucAtt)/10));
        KernelParameters{3}(2:end,1)=RangeValues(1,1,1:length(KernelParameters{3}(2:end,1)));
    end
end
cla(app.OutputPlot,'reset')

%%   single assessment
if OperationIndex==1
    % SSM
    set(app.StatusBar,'Text','Status: Running Assessment SigmaW(L1), SigmaVI, Uncertinty, Speed & Acc. (Two Level)');
    pause(0.1)
    if ~isempty(app.TrueInspectorsData)
        plot(app.OutputPlot,OptBoundsData(2,1)*ones(length(app.TrueInspectorsData{1}),1),app.TrueInspectorsData{1},'o');
        grid(app.OutputPlot,'on');
        ylabel(app.OutputPlot,'Inspector True $\sigma(I_i)$','interpreter','latex')
        xlabel(app.OutputPlot,'Inspector Estimated $\sigma(I_i)$','interpreter','latex')
        drawnow;
    end
    AssessmentNum6();
    set(app.StatusBar,'Text','Status: SSM model parameters estimation is complete');
else
    % SSM-KR
    set(app.StatusBar,'Text','Status: Running Assessment SigmaW(L1), SigmaVI, Uncertinty, Speed & Acc. (Two Level)');
    pause(0.1)
    if ~isempty(app.TrueInspectorsData)
        plot(app.OutputPlot,OptBoundsData(2,1)*ones(length(app.TrueInspectorsData{1}),1),app.TrueInspectorsData{1}(:,1),'o');
        grid(app.OutputPlot,'on');
        ylabel(app.OutputPlot,'Inspector True $\sigma(I_i)$','interpreter','latex')
        xlabel(app.OutputPlot,'Inspector Estimated $\sigma(I_i)$','interpreter','latex')
        drawnow;
    end
    AssessmentNum6();
    set(app.StatusBar,'Text','Status: SSM-KR model parameters estimation is complete');
end
