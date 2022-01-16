function [tot_ll,InterventionMu_Network,InterventionVar_Network]= ...
    opt_intervention_cat(Cindex,InspectorsID,InspectorsData,Re,InspecBiase,y,...
    param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
    InterventionVector,model_i,RegressionModel,AllAtt,InterventionMu_param)
tot_ll=0;
Kernel_l=RegressionModel.Kernel_l;
X_ControlPoints=RegressionModel.X_ControlPoints;
KernelType=RegressionModel.KernelType;
InirilizedEx=RegressionModel.InirilizedEx;
InirilizedVar=RegressionModel.InirilizedVar;
Var_w0=RegressionModel.Sigma_W0^2;

% Prior Knowledge
InterventionMu_Network=InterventionMu_param;
InterventionVar_Network=diag([ModelParamLocal(4)^2 ModelParamLocal(5)^2 ModelParamLocal(6)^2]);%[5^2 0.3^2 0.05^2]);%

Q_r=diag([ModelParamLocal(1)^2 ModelParamLocal(2)^2 ModelParamLocal(3)^2]);%diag([2^2 0.1^2 0.01^2]);
Q_r=blkdiag(Q_r,Q_r);

for i=1:length(Cindex)
    for j=1:length(InspectorsID(1,i,:))
        InspID=find(InspectorsID(1,Cindex(i),j)==InspectorsData{1}(:,1));
        if ~isempty(InspID)
            Re(1,Cindex(i),j)=InspectorsData{1}(InspID,3)^2;
            InspecBiase(1,Cindex(i),j)=InspectorsData{1}(InspID,2)^2;
        end
    end
    % Initial values
    init_x(1)=y(1,Cindex(i),2);
    init_x(3:6,1)=zeros(4,1);
    init_V(1,1)=max(param(3).^2,Re(1,Cindex(i),2));
    init_V(3,3)=param(5).^2;
    init_V(4:6,4:6)=InterventionVar_Network;
    % Kernel Regression
    if model_i==2
        Kr=1;
        for im=1:length(KernelType)
            Krv=KernelFun(AllAtt(i,im),X_ControlPoints(:,im),Kernel_l(im),string(KernelType(im)));
            Kr=Kr.*Krv;
        end
        AKr=Kr./sum(Kr,2);
        init_x(2,:)=AKr*InirilizedEx;
        init_V(2,2)=(AKr*InirilizedVar(:,:,1)*AKr'+Var_w0);
    else
        init_x(2,:)=0;
        [MAxCondition]=RevSpaceTransform(Ncurve,100,100.001,25);
        DifferenceObs=MAxCondition-init_x(1);
        init_V(2,2)=(param(4).^2).*(DifferenceObs)+(param(6).^2);
    end
    
    [~, ~, loglikelihood, ~, ~, InterventionMu_Network,...
        InterventionVar_Network]=HandCoded_KF_Network(y(1,Cindex(i),:),A,F,Q,Q_r,...
        Re(1,Cindex(i),:),InspecBiase(1,Cindex(i),:),param,init_x,init_V,Ncurve,ConstrainedKF,...
        InterventionCheck,InterventionVector(1,Cindex(i),:),...
        InterventionMu_Network,InterventionVar_Network);
    tot_ll= tot_ll + loglikelihood;
end
end