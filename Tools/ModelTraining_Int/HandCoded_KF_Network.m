function [Ex, VarKF, loglikely, Exsmooth, Vsmooth, InterventionMu,...
    InterventionVar,MUSpeed,VarSpeed,MUsSpeed,VarsSpeed,NCDF]=HandCoded_KF_Network(y,A,F,Q,Q_r,R,Bias,param,init_x,init_V,Ncurve,...
    ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu_Net,InterventionVar_Net)
if init_x(2)==0
    % Filter
    [Ex,VarKF] = KalmanFilterFun(y, A, F, Q, Q_r, R, Bias, init_x, init_V,ConstrainedKF,...
                    InterventionCheck,InterventionVector,InterventionMu_Net,InterventionVar_Net,param);
    % Smoother
    TotalTimeSteps=length(Ex(1,:));
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [Exsmooth, ~,~,~]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,Q_r,A,TotalTimeSteps,...
        ConstrainedKF,InterventionCheck,InterventionVector,param);
    InitialCond=Exsmooth(1,2);
    OriginalValues=100;
    MAxCondition=OriginalValues;
    [Mtrv]=RevSpaceTransform(Ncurve,InitialCond,100.0001,25);
    DfferenceObs=MAxCondition(end)-Mtrv;
    % Update Initial State Values
    % Initial variance
    init_V(1,1)=max(R(2),param(2)^2);
    init_V(2,2)=param(4)^2*DfferenceObs+param(6)^2;
    if isnan(init_V(2,2))
        init_V(2,2)
    end
    init_V(3,3)=param(5)^2.;%init_V(2,2)*param(6)^2;
end
% Filter
[Ex,VarKF,~,loglikely,MUSpeed,VarSpeed] = KalmanFilterFun(y, A, F, Q, Q_r, R, Bias, init_x, init_V,...
    ConstrainedKF,InterventionCheck,InterventionVector,...
    InterventionMu_Net,InterventionVar_Net,param);
% Smoother
TotalTimeSteps=length(Ex(1,:));
ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
[Exsmooth, Vsmooth, InterventionMu, InterventionVar, MUsSpeed, VarsSpeed] = KalmanSmootherFun(Ex,...
    VarKF,ExF,VarF,Q,Q_r,A,TotalTimeSteps,ConstrainedKF,InterventionCheck,...
    InterventionVector,param);
IntTime=find(InterventionVector);
for j=0:length(Exsmooth(1,IntTime:end))-1
    NCDF(j+1)=normcdf(Exsmooth(1,IntTime-1),Exsmooth(1,IntTime+j),sqrt(Vsmooth(1,1,IntTime+j)));
end
NPMF=diff(NCDF);
[~,IndvYear]=max(NPMF);
end



function [ExSmooth, VarSmooth,InterventionMu,InterventionSigma,...
    SpeedTracking_mu,SpeedTracking_var]=KalmanSmootherFun(x,Var,ExSmooth,VarSmooth,...
    Q,Q_r, A,TotalTimeSteps,ConstrainedKF,InterventionCheck,InterventionVector,...
    param,varargin)
args = varargin;
if ~isempty(args)
    Ending=args{1};
else
    Ending=1;
end

for i=TotalTimeSteps-1:-1:Ending
    if length(InterventionVector)>1
        if InterventionVector(i+1)
            A(1:3,4:6)=eye(3);
            % added
%             PredVal=A*x(:,i);
%             OCond=RevSpaceTransform(4,PredVal(1));
%             DCond=(100-OCond);
%             Q_r([2,5],[2,5])=Q_r([2,5],[2,5]).*DCond+param^2*diag([1 1]);
            % added
            Q_R=Q_r;
        else
            A(1:3,4:6)=zeros(3);
            Q_R=zeros(6);
        end
    end
    Xpred=A*x(:,i);
    Vpred=A*Var(:,:,i)*A' + Q + Q_R;
    %%
%     if InterventionVector(i+1)
%         if Xpred(2)>ExSmooth(2,i+1)
%             d=[0;-ExSmooth(2,i+1)];
%             D=[0 0 0 0 1 0;0 0 0 0 1 0];
%             [x(:,i),Var(:,:,i)]=KFConstraintsHandlingSeq(x(:,i),Var(:,:,i),D,d,1);
%             Xpred=A*x(:,i);
%             Vpred=A*Var(:,:,i)*A' + Q + Q_R;
%         end
%     end
    %%
    J=Var(:,:,i)*A'*pinv(Vpred);
    ExSmooth(:,i)=x(:,i)+J*(ExSmooth(:,i+1)-Xpred);
    VarSmooth(:,:,i)=Var(:,:,i)+J*(VarSmooth(:,:,i+1)-Vpred)*J';
    ViIndex=find((2*sqrt(VarSmooth(2,2,i))+ExSmooth(2,i))>0);
    if InterventionVector(i+1)
        SpeedTracking_mu(1,1)=ExSmooth(2,i);
        SpeedTracking_mu(2,1)=ExSmooth(5,i);
        SpeedTracking_var(1,1)=VarSmooth(2,2,i);
        SpeedTracking_var(2,1)=VarSmooth(5,5,i);
    end
    KSConstraints();
    if InterventionVector(i+1)
        SpeedTracking_mu(3,1)=ExSmooth(2,i);
        SpeedTracking_mu(4,1)=ExSmooth(5,i);
        SpeedTracking_var(3,1)=VarSmooth(2,2,i);
        SpeedTracking_var(4,1)=VarSmooth(5,5,i);
    end
    if InterventionVector(i+1)
        InterventionMu=ExSmooth(4:6,i);
        InterventionSigma=VarSmooth(4:6,4:6,i);
    end

end
end

