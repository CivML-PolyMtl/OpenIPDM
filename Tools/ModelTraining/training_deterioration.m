
OptBoundsData=app.AllElementsParameters{Index,2}{1};
if model_i==1
    IncludeStructuralAtt=0;
    NumAttributes=[];
else
    IncludeStructuralAtt=1;
    NumAttributes=app.AllElementsParameters{Index,2}{5};
    KernelParameters{1}=app.AllElementsParameters{Index,2}{2}(2:end);
    KernelParameters{2}=app.AllElementsParameters{Index,2}{3};
    KernelParameters{3}=app.AllElementsParameters{Index,2}{4};
end
% Read and Orgnize the dataset
if ~exist([FullPath 'InspectionData_' ElementName '.mat'],'file')
    AppLink.TrainingApp=app;
    AppLink.StatusDropDown.Value=app.AllElementsParameters{Index,2}{7};
    InterventionCond=2;% All: 1, without Interventions: 2, Interventions only: 3
    [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond);
%     InterventionCond=3;% All: 1, without Interventions: 2, Interventions only: 3
%     HandleReadData(AppLink,InterventionCond);
    clear('AppLink');
    InspectorsData=load(FullPathInspector);
    InspectorsData=struct2cell(InspectorsData);
    InspectionData=load(FullPathEx);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
else
    InspectionData=load([FullPath 'InspectionData_' ElementName '.mat']);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
    InspectorsData=load([FullPath 'Inspectors_' ElementName '.mat']);
    InspectorsData=struct2cell(InspectorsData);
end
if ~exist([FullPath 'TrainingData_' ElementName '.mat'],'file')
    TestSet=app.TestsetSpinner.Value/100;
    ParComp=2;
    [ElementData]=OrgnizeData(InspectionData,90E10,[0 0],InspectorsData{1},...
        TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes,TestSet,0);
    save([FullPath 'TrainingData_' ElementName '.mat'],'ElementData');
end
ElementData=load([FullPath 'TrainingData_' ElementName '.mat']);
ElementData=struct2cell(ElementData);
ElementData=ElementData{1};
config_opt.elements = app.AllElementsParameters;
config_opt.indexes = [model_i;elem_j];
save([SavePath '/AutoSave_Config_' ElementName '.mat'],'config_opt');
RunTraining();
AllElementsParameters=app.AllElementsParameters(:,1);
AllElementsParameters{Index,2}=PARAMstore;
AllElementsParameters{Index,3}=InspectorsData{1};
AllElementsParameters{Index,4}=RegressionModelstore;
ElementMetaData=load([FullPath 'MetaData_' ElementName '.mat']);
ElementMetaData=struct2cell(ElementMetaData);
ElementMetaData=ElementMetaData{1};
if model_i==2
    AttOrder=load([FullPath 'RegressAtt_' ElementName '.mat']);
    AttOrder=struct2cell(AttOrder);
    AttOrder=AttOrder{1};
    AttOrder=str2double(cellstr(AttOrder));
    AttOrder=AttOrder(1:NumAttributes)+4;
else
    AttOrder=[];
end
AllElementsParameters{Index,5}=AttOrder;
AllElementsParameters{Index,9}=ElementMetaData;
AllElementsParameters{Index,10}=Nstore;