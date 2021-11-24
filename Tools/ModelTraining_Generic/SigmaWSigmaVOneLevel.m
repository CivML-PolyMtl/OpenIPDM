OptLevel=1;
TruncatedSigma=10^-4;
UpdatedInspectorsData{1}=[UpdatedInspectorsData{1}(:,1) ...
    OptBoundsData(2,1)*ones(length(UpdatedInspectorsData{1}),1)];

UpdatedInspectorsData{1}(:,3)=OptBoundsData(2,1);
UpdatedInspectorsData{1}(:,2) = 0;

% OperationIndex=find(strcmp(app.SelectOptAlgorithm.Value,...
%     app.SelectOptAlgorithm.Items));

if OptimizationAlgorithmIndex==1 || OptimizationAlgorithmIndex==2
    [Qparam,~,~,fx_NR]=Newton_Raphson_par(@(PARAM) AnalysisObjective(MdataEngy,...
        UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
        x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
        DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
        GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,GradCompute)...
        ,PARAM,'log_transform','no','output','original','laplace','no',...
        'convergence_tol',1E-3,'bounds',boundsQ);
%     Qparam=PARAM;fx_NR=10^-10;
    UpdatedInspectorsData{1}=[UpdatedInspectorsData{1}(:,1)...
        OptBoundsData(2,1)*ones(length(UpdatedInspectorsData{1}),1)...
        Qparam(2)*ones(length(UpdatedInspectorsData{1}),1)];
    LogLik=fx_NR(end); 
end
