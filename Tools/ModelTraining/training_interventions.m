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
if ~exist([FullPath 'InspectionData_Intervention_' ElementName '.mat'],'file')
    AppLink.TrainingApp=app;
    AppLink.StatusDropDown.Value=app.AllElementsParameters{Index,2}{7};
%     InterventionCond=2;% All: 1, without Interventions: 2, Interventions only: 3
%     [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond);
    InterventionCond=3;% All: 1, without Interventions: 2, Interventions only: 3
    [FullPathEx,FullPathInspector]=HandleReadData(AppLink,InterventionCond);
    clear('AppLink');
%     InspectorsData=load(FullPathInspector);
%     InspectorsData=struct2cell(InspectorsData);
    InspectionData=load(FullPathEx);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
else
    InspectionData=load([FullPath 'InspectionData_Intervention_' ElementName '.mat']);
    InspectionData=struct2cell(InspectionData);
    InspectionData=InspectionData{1};
%     InspectorsData=load([FullPath 'Inspectors_' ElementName '.mat']);
%     InspectorsData=struct2cell(InspectorsData);
end
if ~exist([FullPath 'TrainingData_Intervention_' ElementName '.mat'],'file')
    TestSet=0;
    ParComp=1;
    [ElementData]=OrgnizeData(InspectionData,90E10,[0 0],InspectorsData{1},...
        TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes,TestSet,1);
    save([FullPath 'TrainingData_Intervention_' ElementName '.mat'],'ElementData');
else
    ElementData=load([FullPath 'TrainingData_Intervention_' ElementName '.mat']);
    ElementData=struct2cell(ElementData);
    ElementData=ElementData{1};
end



Yearly=ElementData.yearlyS;
InspectorIDLabel=ElementData.InspectorLabelS;


%% KR
if model_i==2
    AttStruc=reshape((ElementData.StrucAtt),size(ElementData.StrucAtt,2),size(ElementData.StrucAtt,3));
    AvgObs=reshape((ElementData.init_x),size(ElementData.init_x,2),1);
    AllAtt=[AttStruc AvgObs];
else
    RegressionModel = [];
    AllAtt = [];
end

param = Qparam;

% time step
dt=1;

% transition for the kinimatic model
Ax1=[1 dt dt^2/2;
    0 1 dt;
    0 0 1];

A=blkdiag(Ax1,eye(3));

% inspection data
y = (ElementData.YS);

% Observation matrix
F =[1 zeros(1,5)];

% Noise
Re = ElementData.ReS;
Be = ElementData.InpecBiaseS;

Q_ki = param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;
    (dt^4)/8 (dt^3)/3 (dt^2)/2;
    (dt^3)/6 (dt^2)/2 dt];
Q=blkdiag(Q_ki,zeros(3));

% Apply Constraints
ConstrainedKF = 1;

% Intervention Anlayese
InterventionCheck = 1;

% Intervention Vector
InterventionVector = (ElementData.InterventionVector);

init_x=zeros(6,1);
[~,MAxCondition]=SpaceTransformation(Ncurve,100,100,25);

IntervType = (ElementData.InterventionType);
IntervType = fix(IntervType./10.^fix(log10(IntervType)));

OptBounds_int = {OptBoundsInt(1,:) OptBoundsInt(2,:) OptBoundsInt(3,:) ...
    OptBoundsInt(4,:) OptBoundsInt(5,:) OptBoundsInt(6,:)};

for IndexTypeVal=1:3
    Cindex = find(IntervType==IndexTypeVal);
    if ~isempty(Cindex)
        InspectorsID=(ElementData.InspectorLabelS);
        [ModelParamLocal,~,~,fx_NR]=Newton_Raphson(@(ModelParamLocal)...
            opt_intervention(Cindex,InspectorsID,InspectorsData,Re,Be,y,...
            param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
            InterventionVector,model_i,RegressionModel,AllAtt)...
            ,ModelParamLocal,'log_transform','no','output','original','laplace','no',...
            'convergence_tol',1E-3,'bounds',OptBounds_int);
        [~,InterventionMu_Network,InterventionVar_Network]= ...
            opt_intervention(Cindex,InspectorsID,InspectorsData,Re,Be,y,...
            param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
            InterventionVector,model_i,RegressionModel,AllAtt);
        int_Ex{1,IndexTypeVal} = InterventionMu_Network;
        int_Ex{2,IndexTypeVal} = InterventionVar_Network;
        int_param {1,IndexTypeVal} = ModelParamLocal;
    else
        int_Ex{1,IndexTypeVal} = nan(3,1);
        int_Ex{2,IndexTypeVal} = nan(3);
        int_param {1,IndexTypeVal} = nan(1,6);
    end
end