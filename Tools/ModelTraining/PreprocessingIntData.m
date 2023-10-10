OptBoundsInt = [];
if size(app.AllElementsParameters,2)<2
    OptBoundsData = [];
elseif ~isempty(find(~cellfun(@isempty,app.AllElementsParameters(:,2))))
    [~,ind_optbounds] = max(cellfun(@length,app.AllElementsParameters(:,2)));
    app.AllElementsParameters{Index,2}{8} = app.AllElementsParameters{ind_optbounds, 2}{8};
    OptBoundsInt = app.AllElementsParameters{ind_optbounds, 2}{8}(:,2:3);
end
% if the user already defines the bounds 
if size(app.AllElementsParameters,2)>1 && ~isempty(app.AllElementsParameters{Index,2})
    OptBoundsInt = app.AllElementsParameters{Index,2}{8}(:,2:3);
end
if isempty(OptBoundsInt)
    load(sprintf(['%s/Tools/ModelTraining/InitialParameterValues/SoftwareInitialInt_net.mat'],pwd));
    app.AllElementsParameters{Index,2}{8} = SoftwareInitialInt;
    OptBoundsInt=app.AllElementsParameters{Index,2}{8}(:,2:3);
end

if model_i==1
    IncludeStructuralAtt=0;
else
    IncludeStructuralAtt=1;
    NumAttributes = app.AllElementsParameters{Index,2}{5}; 
    % TAGI parameters
    if isempty(AnnModel)
        AnnModel.batch_size = 2;%TAGI_params{2};
        AnnModel.max_epoches = 30;%TAGI_params{1};
        AnnModel.cont_training = 0;
        AnnModel.evaluate_mode = 1;
        AnnModel.theta_nn = [];
    end
end

% Read and Orgnize the dataset
if ~exist([FullPath 'InspectionData_Intervention_' erase(ElementName,"/") '.mat'],'file')
    AppLink.TrainingApp=app;
    AppLink.StatusDropDown.Value=app.AllElementsParameters{Index,2}{7};
%     InterventionCond=2;% All: 1, without Interventions: 2, Interventions only: 3
%     [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond);
    InterventionCond=3;% All: 1, without Interventions: 2, Interventions only: 3
    [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond,IncludeStructuralAtt);
    clear('AppLink');
%     InspectorsData=load(FullPathInspector);
%     InspectorsData=struct2cell(InspectorsData);
    if ~isempty(FullPathEx)
        InspectionData=load(FullPathEx);
        InspectionData=struct2cell(InspectionData);
        InspectionData=InspectionData{1};
    else
        InspectionData = [];
    end
else
    InspectionData=load([FullPath 'InspectionData_Intervention_' erase(ElementName,"/") '.mat']);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
%     InspectorsData=load([FullPath 'Inspectors_int_' ElementName '.mat']);
%     InspectorsData=struct2cell(InspectorsData);
end
if ~exist([FullPath 'TrainingData_Intervention_' erase(ElementName,"/") '.mat'],'file')
    TestSet=0;
    ParComp=1;
    if strcmp(app.Switch.Value, 'Real')
        syn_flag = 0;
    else
        syn_flag = 1;
    end
    if ~isempty(InspectionData)
        [ElementData, Success]=OrgnizeData(InspectionData,90E10,[0 0],InspectorsData{1},...
            TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes,TestSet,1,syn_flag);
    else
        Success = 0;
    end
    if Success
        save([FullPath 'TrainingData_Intervention_' erase(ElementName,"/") '.mat'],'ElementData');
    end
else
    ElementData=load([FullPath 'TrainingData_Intervention_' erase(ElementName,"/") '.mat']);
    ElementData=struct2cell(ElementData);
    ElementData=ElementData{1};
    Success = 1;
end