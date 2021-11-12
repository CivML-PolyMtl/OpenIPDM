%% SSM
dt=1;                                                                       % time step
Ak=[1 dt dt^2/2;0 1 dt;0 0 1];                                              % transition for the kinimatic model
A=blkdiag(Ak,diag(ones(3,1)));                                                       % transition matrix (with intervention)
C =[1 zeros(1,5)];                                                          % Observation matrix
Q_ki = PriorParam(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;
                    (dt^4)/8 (dt^3)/3 (dt^2)/2;
                    (dt^3)/6 (dt^2)/2 dt];                                  % Process error (covariance)
Q=blkdiag(Q_ki,zeros(3));                                                   % Process error (full matrix)                              
for j=1:length(Inspectors(:))                                               % Observation error (variance)
    InspID=find(Inspectors(j)==InspectorsParam(:,1));
    if ~isempty(InspID)
        Re(1,j)=InspectorsParam(InspID,3)^2;
        if app.HideInspectorIDMenu.Checked
            InspectorIDLabel_y(j)=InspID;
        else
            InspectorIDLabel_y(j)=InspectorsParam(InspID,1);
        end
    else
        Re(1,j)=mean(InspectorsParam(:,3))^2;
        InspectorIDLabel_y(j)=0;
    end
end
R_values=nan(TimeWindow,1);
R_values(iaY)=Re;
Re=R_values;
ConstrainedKF = 1;                                                          % Apply PDF Truncation

%% Interventions
IntervType=round(InterventionType/1000);                                    % Intervention Type
if InterventionType~=0
    for InType=1:length(IntervType)
        R_param=IntParam{1,IntervType}; % col: InterventionType
        Q_r1=diag([R_param(1)^2 R_param(2)^2 R_param(3)^2]);                         % Local error (covariance)
        Q_r{InType}=blkdiag(Q_r1,Q_r1);  
    end
else
    R_param=nan(1,6);
    Q_r1=diag([R_param(1)^2 R_param(2)^2 R_param(3)^2]);                         % Local error (covariance)
    Q_r{1}=blkdiag(Q_r1,Q_r1);  
end
for InType=1:length(IntervType)
    if isnan(Q_r{InType}(1,1)) && IntervType(InType)~=0
        R_param=NAParam{ImportanceInd,IntervType}; % col: InterventionType
        Q_r1=diag([R_param(1)^2 R_param(2)^2 R_param(3)^2]);                         % Local error (covariance)
        Q_r{InType}=blkdiag(Q_r1,Q_r1);
    end
end
FindType=0;
InterventionCheck = [];
InterventionVector = zeros(TimeWindow,1);
IintY=find(sum(InterventionYear==YearSteps,1));
InterventionVector(IintY)=1;
                                                     % Local error (full matrix)

if max(IntervType)~=0 && sum(InterventionVector)>0 && min(find(InterventionVector))>1
    for InType=1:length(IntervType)
        InterventionVar_Network{InType}=VarIntParam{IntervType};                         % Prior effect (covariance)                            
        InterventionMu_Network{InType}=ExIntParam{IntervType};
        % Default Int value
        if sum(isnan(InterventionVar_Network{InType}))>0
            InterventionVar_Network{InType}=NAIntParam{ImportanceInd+2,IntervType};                                       % Prior effect (covariance)                            
            InterventionMu_Network{InType}=NAIntParam{ImportanceInd,IntervType};
        end
    end
else
    InterventionVar_Network{1}=diag(ones(3,1));                                       % Prior effect (covariance)                            
    InterventionMu_Network{1}=nan(3,1);
end

%% AutoCorrect Check
if strcmp(app.AutocorrectElem.Value,'On')
    y_Data_before(:,1)=YearTotal;
    y_Data_before(:,2)=y_Data;
    [YearTotal,y_Data,Re,IntVectorCorrect,InspectorIDLabel_y,InterventionYear]=...
        da_AutoCorrect(app,YearTotal,y_Data,Re,InterventionVector,...
        InspectorIDLabel_y,InterventionDuration,BudgetDuration,InterventionYear);
    if sum(IntVectorCorrect)&& ~sum(InterventionVector)
        InterventionVector = zeros(length(YearTotal),1);
        IintY=find(InterventionYear==YearTotal);
        InterventionVector(IintY)=1;
        FindType=1;
    elseif sum(InterventionVector)
        InterventionVector = zeros(length(YearTotal),1);
        [~,IintY]=intersect(YearTotal,InterventionYear);
        InterventionVector(IintY)=1;
    end
else
    y_Data_before=[];
end

%% Initial state
init_x=zeros(6,1);
init_V=eye(6);
init_x(1)=y_Data(2);
init_V(1,1)=max(PriorParam(3).^2,Re(2));
init_V(3,3)=PriorParam(5).^2;
init_V(4:6,4:6)=InterventionVar_Network{1};


%% KR
if ~isempty(RegressionModel)
    Kernel_l=RegressionModel.Kernel_l;
    X_ControlPoints=RegressionModel.X_ControlPoints;
    KernelType=RegressionModel.KernelType;
    InirilizedEx=RegressionModel.InirilizedEx;
    InirilizedVar=RegressionModel.InirilizedVar;
    Var_w0=RegressionModel.Sigma_W0^2;
    AttStruc=StructuralAttributes(:,SA_Index-4);
    IndNaNCheck=find(isnan(AttStruc));
    if ~isempty(IndNaNCheck)
        AttStruc(IndNaNCheck)=mean(X_ControlPoints(:,IndNaNCheck));
    end
    if sum(find(SA_Index==14))>0
        Lind=find(SA_Index==14);
        ZLind=find(AttStruc(:,Lind)==0);
        AttStruc(ZLind,Lind)=1;
        AttStruc(:,Lind)=log(AttStruc(:,Lind));
    end
    y_Data_ind=find(~isnan(y_Data));
    if length(y_Data_ind)>2
        AvgObs=mean(y_Data(y_Data_ind(1:3)));
    else
        AvgObs=mean(y_Data(y_Data_ind));
    end
    AllAtt=[AttStruc AvgObs];
    
    Kr=1;
    for im=1:length(KernelType)
        Krv=KernelFun(AllAtt(im),X_ControlPoints(:,im),Kernel_l(im),...
            string(KernelType(im)));
        Kr=Kr.*Krv;
    end
    AKr=Kr./sum(Kr,2);
    init_x(2,:)=AKr*InirilizedEx;
    init_V(2,2)=(AKr*InirilizedVar(:,:,1)*AKr'+Var_w0);
    if init_x(2,:)>0 
        D=[1;1];
        d=[-5;0];
        [init_x(2,1),init_V(2,2,1)]=KFConstraintsHandlingSeq(init_x(2,1),init_V(2,2,1),D,d,1);
    end
else
    init_x(2,:)=0;
    [Mtrv]=RevSpaceTransform(Ncurve,init_x(1),100,25);
    DfferenceObs=100-Mtrv;
    init_V(2,2)=PriorParam(4)^2*DfferenceObs+PriorParam(6)^2;
end



