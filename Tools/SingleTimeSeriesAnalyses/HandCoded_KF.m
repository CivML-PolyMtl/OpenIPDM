function [Ex, VarKF, param, loglikely, Exsmooth, Vsmooth]=HandCoded_KF(y,A,...
    F,Q,R,param,init_x,init_V,Ncurve,ConstrainedKF,KFCase,InterventionCheck,...
    InterventionVector,InterventionMu,InterventionSigma)
if KFCase==1
[Ex,VarKF] = KalmanFilterFun(y, A, F, Q, R, init_x, init_V,ConstrainedKF,...
    InterventionCheck,InterventionVector,InterventionMu,InterventionSigma,...
    param);
TotalTimeSteps=length(Ex(1,:));
if InterventionCheck
    TimeBeforeIntervention=find(InterventionVector);
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [ExsmoothPI, VsmoothPI]=KalmanSmootherFun(Ex,VarKF, ExF,VarF,Q,A,TotalTimeSteps,...
    ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
    InterventionSigma,param,KFCase,TimeBeforeIntervention);
    ExF(:,TimeBeforeIntervention-1)=Ex(:,TimeBeforeIntervention-1);
    VarF(:,:,TimeBeforeIntervention-1)=VarKF(:,:,TimeBeforeIntervention-1);
    [Exsmooth, Vsmooth]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,A,TimeBeforeIntervention-1,...
        ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
        InterventionSigma,param,KFCase);
    Exsmooth(:,TimeBeforeIntervention:end)=ExsmoothPI(:,TimeBeforeIntervention:end);
    Vsmooth(:,:,TimeBeforeIntervention:end)=VsmoothPI(:,:,TimeBeforeIntervention:end);
else
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [Exsmooth, Vsmooth]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,A,TotalTimeSteps,...
        ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
        InterventionSigma,param,KFCase);
end

InitialCond=Exsmooth(1,2);
OriginalValues=100;
MAxCondition=OriginalValues;
[Mtrv]=RevSpaceTransform(Ncurve,InitialCond);
DfferenceObs=MAxCondition(end)-Mtrv;

% Update Initial State Values
% Initial expected value
init_x(1)=y(2);%nanmean(y);%
init_x(2)=param(5)*DfferenceObs;
init_x(3)=param(6);%*init_x(2);
% Initial variance
init_V(1,1)=max(R(2),param(2)^2);
init_V(2,2)=param(3)^2*DfferenceObs+param(7)^2;
init_V(3,3)=param(4)^2.;%init_V(2,2)*param(6)^2;


[Ex,VarKF,~,loglikely] = KalmanFilterFun(y, A, F, Q, R, init_x, init_V,...
    ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,InterventionSigma,param);
if InterventionCheck
    TimeBeforeIntervention=find(InterventionVector);
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [ExsmoothPI, VsmoothPI]=KalmanSmootherFun(Ex,VarKF, ExF,VarF,Q,A,TotalTimeSteps,...
    ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
    InterventionSigma,param,KFCase,TimeBeforeIntervention);
    ExF(:,TimeBeforeIntervention-1)=Ex(:,TimeBeforeIntervention-1);
    VarF(:,:,TimeBeforeIntervention-1)=VarKF(:,:,TimeBeforeIntervention-1);
    [Exsmooth, Vsmooth]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,A,TimeBeforeIntervention-1,...
        ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
        InterventionSigma,param,KFCase);
    Exsmooth(:,TimeBeforeIntervention:end)=ExsmoothPI(:,TimeBeforeIntervention:end);
    Vsmooth(:,:,TimeBeforeIntervention:end)=VsmoothPI(:,:,TimeBeforeIntervention:end);
else
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [Exsmooth, Vsmooth]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,A,TotalTimeSteps,...
        ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
        InterventionSigma,param,KFCase);
end
% ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
% VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
% [Exsmooth, Vsmooth]=KalmanSmootherFun(Ex,VarKF, ExF,VarF,Q,A,TotalTimeSteps,...
%     ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu,...
%     InterventionSigma,param,TimeBeforeIntervention);
elseif KFCase==2
    yind=find(~isnan(y));
    init_x(1)=nanmean(y(yind(1:3)));%
    init_x(2)=param(5);
    init_x(3)=param(6);%*init_x(2);
    % Initial variance
    init_V(1,1)=max(R(2),param(2)^2);
    init_V(2,2)=param(3)^2;
    init_V(3,3)=param(4)^2;%init_V(2,2)*param(6)^2;
    % Equations
%     TransitionFunction=@(Ex) A*Ex+[-1./(1+exp(Ex(2,:)))
%         zeros(1,length(Ex(2,:)));
%         zeros(1,length(Ex(2,:)))];
    dt=1;
    A=eye(3);
    npDotX=4;
    npDDotX=4;
    TransitionFunction=@(Ex) A*Ex+[TransFun(npDotX,Ex(2,:),0,-2)+TransFun(npDotX,Ex(3,:),0,-2)*dt^2/2; 
        TransFun(npDDotX,Ex(3,:),0,-2)*dt;    
        zeros(1,6)];
    MeasurementFunction=@(Ex) F*Ex;
    
    [Ex,VarKF,Std_KF,loglikely]=CubtureKalmanFilterFun(TransitionFunction,...
        MeasurementFunction, y, Q, R, init_x', init_V);
%     [Ex,VarKF,Std_KF,loglikely]=CubtureKalmanFilter(TransitionFunction,...
%         MeasurementFunction, y, Q, R, F, init_x', init_V);
    [Exsmooth, Vsmooth]=CubatureKalmanSmootherFun(TransitionFunction, Ex,...
        VarKF, Std_KF, Q);
end
end

function [ExpectedCond, Variance, SigmaForecast, loglik]=KalmanFilterFun(y,...
    A, F, Q, R, init_x, init_V,ConstrainedKF,InterventionCheck,...
    InterventionVector,InterventionMu,InterventionSigma,param)
T=length(y);

% Equations
UpdateExpected=@(EX,EY,sigma_XY,sigma_obs,y) EX+sigma_XY*sigma_obs^-1*(y-EY);
UpdateVariance=@(sigma_X2,sigma_obs,sigma_XY) sigma_X2-sigma_XY*sigma_obs^-1*sigma_XY';

ExpectedCond(:,1)=init_x;
Variance(:,:,1)=init_V;
SigmaForecast(:,:,1)=sqrt(diag(init_V));
loglik=zeros(T,1);
for t=0:T-1
    if t==0
        UPvariance=Variance;
        S=F*UPvariance*F'+R(t+1);
        sigmaXY=UPvariance*F';
    else
        if length(InterventionVector)>1
            if InterventionVector(t+1) 
                A(1:3,4:6)=eye(3);
                ExpectedCond(4:6,t)=InterventionMu;
            else
                A(1:3,4:6)=zeros(3);
            end
        end
        UPvariance=A*Variance(:,:,t)*A'+Q;
        S=F*UPvariance*F'+R(t+1);
        sigmaXY=UPvariance*F';
    end
    if t~=0
        if ~isnan(y(t+1))
        ExpectedCond(:,t+1)=UpdateExpected(A*ExpectedCond(:,t),F*(A*ExpectedCond(:,t)),sigmaXY,S,y(t+1)); 
        Variance(:,:,t+1)=UpdateVariance(UPvariance,S,sigmaXY);
        else
            ExpectedCond(:,t+1)=A*ExpectedCond(:,t);
            Variance(:,:,t+1)=UPvariance;
        end
        if length(InterventionVector)>1
            if InterventionVector(t+1) 
                %ExpectedCond(2:3,t+1)=0;
                [Mtrv]=RevSpaceTransform(3,ExpectedCond(1,t+1));
                DfferenceObs=100-Mtrv;
                VarSpeed=param(3)^2*DfferenceObs+param(7)^2;
                C=[0 1 zeros(1,4);0 0 1 zeros(1,3)];
                S_speed=C*UPvariance*C'+diag([0.01, 0.001]);
                sigmaXY_speed=UPvariance*C';
                A(1:3,4:6)=zeros(3);
                ExpectedCond(:,t+1)=UpdateExpected(A*ExpectedCond(:,t+1),C*(A*ExpectedCond(:,t+1)),sigmaXY_speed,S_speed,[-0.2,0]');
                Variance(:,:,t+1)=UpdateVariance(UPvariance,S_speed,sigmaXY_speed);
                
            end
        end
        if ConstrainedKF && ExpectedCond(2,t+1)+2*sqrt(Variance(2,2,t+1))>0
            d=[-50;0];
            if InterventionCheck
                D=[0 1 0 0 0 0;0 1 0 0 0 0];
            else
                D=[0 1 0;0 1 0];
            end
            [ExpectedCond(:,t+1),Variance(:,:,t+1)]=KFConstraintsHandlingSeq(ExpectedCond(:,t+1),Variance(:,:,t+1),D,d,1);
        end
        SigmaForecast(:,t+1)=sqrt(diag(Variance(:,:,t+1)));
        if ~isnan(y(t+1))
            Err=y(t+1)-F*(A*ExpectedCond(:,t+1));
            loglik(t+1)=log(normpdf(Err,0,sqrt(S)));
        else
            loglik(t+1)=NaN;
        end
    end
end
loglik=sum(loglik(~isnan(loglik)));
end 

function [ExSmooth, VarSmooth]=KalmanSmootherFun(x,Var,ExSmooth,VarSmooth,...
    Q,A,TotalTimeSteps,ConstrainedKF,InterventionCheck,InterventionVector,...
    InterventionMu,InterventionSigma,param,KFCase,varargin)
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
        else
            A(1:3,4:6)=zeros(3);
        end
    end
    Xpred=A*x(:,i);
    Vpred=A*Var(:,:,i)*A' + Q;
    J=Var(:,:,i)*A'*inv(Vpred);
    ExSmooth(:,i)=x(:,i)+J*(ExSmooth(:,i+1)-Xpred);
    VarSmooth(:,:,i)=Var(:,:,i)+J*(VarSmooth(:,:,i+1)-Vpred)*J';
    ViIndex=find((2*sqrt(VarSmooth(2,2,i))+ExSmooth(2,i))>0);
    if ~isempty(ViIndex) && ConstrainedKF && KFCase==1
        d=[-50;0];
        if InterventionCheck
            D=[0 1 0 0 0 0;0 1 0 0 0 0];
        else
            D=[0 1 0;0 1 0];
        end
        [ExSmooth(:,i),VarSmooth(:,:,i)]=KFConstraintsHandlingSeq(ExSmooth(:,i),VarSmooth(:,:,i),D,d,1);
    end
end
end