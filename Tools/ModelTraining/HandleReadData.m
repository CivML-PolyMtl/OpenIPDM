function [FullPathEx,FullPath_Inspectors]=HandleReadData(AppLink,DataType,IncludeStructuralAtt)
    SE=AppLink.TrainingApp.SeDatabaseStore;
    ST=AppLink.TrainingApp.StDatabaseStore;
    IP=AppLink.TrainingApp.InspectorsDatabaseStore;
    IN=AppLink.TrainingApp.IntDatabaseStore;
    IB=AppLink.TrainingApp.IntBudgetDataStore;
    ED=AppLink.TrainingApp.ElementDetailsDataStore;
    SelectedElementCat=AppLink.TrainingApp.Tree.SelectedNodes.Text;
    StatusValue=AppLink.StatusDropDown.Value;
    OriginPWD=pwd;
    InterventionCond=DataType;% All: 1, without Interventions: 2, Interventions only: 3
    addpath(sprintf('%s/Tools/ModelTraining/',pwd));
    load(sprintf(['%s/Parameters/ModelTraining_config_elem.mat'],pwd));
    check_exist = find(strcmp(SelectedElementCat,config_elem(:,1)));
    if ~isempty(check_exist)
        if strcmp(config_elem{check_exist,3}, "Material")
            categorical_data_ind = 7; % material column
        elseif strcmp(config_elem{check_exist,3}, "Type")
            categorical_data_ind = 8; % element type column
        else
            categorical_data_ind = 7; % no categorical data
        end
    end
    ExcuteReadData();
end