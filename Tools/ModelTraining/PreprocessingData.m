AnnModel=[];
NumAttributes=0;
% if the user didn't define the bounds at all
if size(app.AllElementsParameters,2)<2
    OptBoundsData = [];
elseif ~isempty(find(~cellfun(@isempty,app.AllElementsParameters(:,2))))
    ind_optbounds = min(find(~cellfun(@isempty,app.AllElementsParameters(:,2))));
    OptBoundsData = app.AllElementsParameters{ind_optbounds, 2}{1};
end
% if the user already defines the bounds 
if size(app.AllElementsParameters,2)>1 && ~isempty(app.AllElementsParameters{Index,2})
    OptBoundsData=app.AllElementsParameters{Index,2}{1};
end
if isempty(OptBoundsData)
    load(sprintf(['%s/Parameters/ModelTraining_SoftwareInitial_net.mat'],pwd));
    load(sprintf(['%s/Parameters/ModelTraining_SoftwareInitialInt_net.mat'],pwd));
    app.AllElementsParameters{Index,2} = [];
    app.AllElementsParameters{Index,2}{1} = BoundsData;
    app.AllElementsParameters{Index,2}{7} = 'Actif Min.';
    app.AllElementsParameters{Index,2}{5} = 0;
    OptBoundsData=app.AllElementsParameters{Index,2}{1};
    app.AllElementsParameters{Index,2}{8} = SoftwareInitialInt;
end

if isempty(app.AllElementsParameters{Index,2})
    app.AllElementsParameters{Index,2}{1} = OptBoundsData;
    app.AllElementsParameters{Index,2}{7} = 'Actif Min.';
    app.AllElementsParameters{Index,2}{5} = NumAttributes;
end

if model_i==1
    IncludeStructuralAtt=0;
else
    IncludeStructuralAtt=1;
    if strcmp(app.Switch.Value, 'Synthetic')
        NumAttributes = 1;
    elseif sum(strcmp(ElementName,{'Chasse-roue / trottoir','Tirants','Revêtement de mur','Voûte / Dalle','Toiture'})) == 1
        NumAttributes = 10;
    elseif sum(strcmp(ElementName, {'Dessous de la dalle/voûte','Mur de tête','Mur en aile','Murs/naiss.voûte/coins infér.'}))== 1
        NumAttributes = 9;
    elseif strcmp(ElementName,'Mur')
        NumAttributes = 6;
    else 
        NumAttributes = 11;
    end
    app.AllElementsParameters{Index,2}{5} = NumAttributes;
    % TAGI parameters
    StructuralAttributes=0;
    if isempty(AnnModel)
        AnnModel.batch_size = 2;%TAGI_params{2};
        AnnModel.max_epoches = 30;%TAGI_params{1};
        AnnModel.cont_training = 0;
        AnnModel.evaluate_mode = 1;
        AnnModel.theta_nn = [];
    end
end
% Read and Orgnize the dataset
if ~exist([FullPath 'InspectionData_' erase(ElementName,"/") '.mat'],'file') && strcmp(app.Switch.Value, 'Real')
    AppLink.TrainingApp=app;
    AppLink.StatusDropDown.Value=app.AllElementsParameters{Index,2}{7};
    InterventionCond=2;% All: 1, without Interventions: 2, Interventions only: 3
    [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond,IncludeStructuralAtt);
    %     InterventionCond=3;% All: 1, without Interventions: 2, Interventions only: 3
    %     HandleReadData(AppLink,InterventionCond);
    clear('AppLink');
    InspectorsData=load(FullPathInspector);
    InspectorsData=struct2cell(InspectorsData);
    InspectionData=load(FullPathEx);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
else
    InspectionData=load([FullPath 'InspectionData_' erase(ElementName,"/") '.mat']);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
    InspectorsData=load([FullPath 'Inspectors_' erase(ElementName,"/") '.mat']);
    InspectorsData=struct2cell(InspectorsData);
end

if ~exist([FullPath 'TrainingData_' erase(ElementName,"/") '.mat'],'file')
    TestSet=app.TestsetSpinner.Value/100;
    ParComp=2;
    if strcmp(app.Switch.Value, 'Real')
        syn_flag = 0;
    else
        syn_flag = 1;
    end
    [ElementData, Success]=OrgnizeData(InspectionData,90E10,[0 0],InspectorsData{1},...
        TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes,TestSet,0,syn_flag);
    if Success
        if IncludeStructuralAtt && ~syn_flag
            ElementData.StrucAtt = cat(3, ElementData.StrucAtt, ElementData.init_x);
            ElementData.ModelValid.StrucAtt = cat(3, ElementData.ModelValid.StrucAtt, ElementData.ModelValid.init_x);
            ElementData.ModelTest.StrucAtt = cat(3, ElementData.ModelTest.StrucAtt, ElementData.ModelTest.init_x);
        end
        save([FullPath 'TrainingData_' erase(ElementName,"/") '.mat'],'ElementData', '-v7.3');
    end
else
    Success = 1;
end