function [FullPathEx,FullPath_Inspectors]=HandleReadData(AppLink,DataType)
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
    ExcuteReadData();
end