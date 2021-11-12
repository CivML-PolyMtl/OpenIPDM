if isempty(app.RegressionModel)
    IncludeStructuralAtt=0;
    NumAttributes=[];
    model_i = 1;
else
    IncludeStructuralAtt=1;
    if SynDataChk
        NumAttributes=length(app.RegressionModel.Kernel_l);
    else
        NumAttributes=length(app.RegressionModel.Kernel_l)-1;
    end
    model_i = 2;
end



TestSet=0;
ParComp=1;
ThisYear=date;
ThisYear=str2double(ThisYear(end-3:end));
TrainingWindow=ThisYear-app.StartYearSpinner.Value;

[ElementData]=OrgnizeDataInt(app.InspectionData,90E10,[0 0],app.InspectorsParams{1}(:,[1 3]),SynDataChk,Ncurve,...
    [],TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes); 
% [ElementData]=OrgnizeData_int(app.InspectionData,90E10,[0 0],app.InspectorsParams{1}(:,[1 3]),...
%     TrainingWindow,ParComp,IncludeStructuralAtt,NumAttributes,TestSet,1);

Yearly=ElementData.yearlyS;
InspectorIDLabel=ElementData.InspectorLabelS;


%% KR
AttStruc=reshape((ElementData.StrucAtt),size(ElementData.StrucAtt,2),size(ElementData.StrucAtt,3));
AvgObs=reshape((ElementData.init_x),size(ElementData.init_x,2),1);
if ~SynDataChk
    AllAtt = [AttStruc AvgObs];
else
    AllAtt = AttStruc;
end

param = app.SSMParams;

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
Re = (ElementData.ReS);
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


OptBounds_int = {cell2mat(OptBoundsInt(1,2:3)) cell2mat(OptBoundsInt(2,2:3)) cell2mat(OptBoundsInt(3,2:3)) ...
    cell2mat(OptBoundsInt(4,2:3)) cell2mat(OptBoundsInt(5,2:3)) cell2mat(OptBoundsInt(6,2:3))};

ModelParamLocal = cell2mat(OptBoundsInt(1:6,1))';

Mu_net = cell2mat(OptBoundsInt(7:end,1))';
Cindex = find(IntervType==str2double(IntCat)/1000);
if ~isempty(Cindex)
    InspectorsID=(ElementData.InspectorLabelS);
    [ModelParamLocal,~,~,fx_NR]=Newton_Raphson_par(@(ModelParamLocal)...
        opt_intervention_cat(Cindex,InspectorsID,app.InspectorsParams,Re,y,...
        param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
        InterventionVector,model_i,app.RegressionModel,AllAtt,Mu_net)...
        ,ModelParamLocal,'log_transform','no','output','original','laplace','no',...
        'convergence_tol',1E-3,'bounds',OptBounds_int);
    [~,InterventionMu_Network,InterventionVar_Network]= ...
        opt_intervention_cat(Cindex,InspectorsID,app.InspectorsParams,Re,y,...
        param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
        InterventionVector,model_i,app.RegressionModel,AllAtt,Mu_net);
    app.int_Ex{1,1} = InterventionMu_Network;
    app.int_Ex{2,1} = InterventionVar_Network;
    app.int_param {1,1} = ModelParamLocal;
else
    app.int_Ex{1,1} = nan(3,1);
    app.int_Ex{2,1} = nan(3);
    app.int_param {1,1} = nan(1,6);
end
