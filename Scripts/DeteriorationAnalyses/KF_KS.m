function [Ex, VarKF, loglikely, Exsmooth, Vsmooth, InterventionMu,...
    InterventionVar,MUSpeed,VarSpeed,MUsSpeed,VarsSpeed,NCDF]=KF_KS(y,A,F,Q,Q_r,R,param,init_x,init_V,Ncurve,...
    ConstrainedKF,InterventionCheck,InterventionVector,InterventionMu_Net,InterventionVar_Net)
if init_x(2)==0
    % Filter
    [Ex,VarKF] = KalmanFilterFun(y, A, F, Q, Q_r, R, init_x, init_V,ConstrainedKF,...
                    InterventionCheck,InterventionVector,...
                    InterventionMu_Net,InterventionVar_Net,param);
    % Smoother
    TotalTimeSteps=length(Ex(1,:));
    ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
    VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
    [Exsmooth, ~,~,~]=KalmanSmootherFun(Ex,VarKF,ExF,VarF,Q,Q_r,A,TotalTimeSteps,...
        ConstrainedKF,InterventionCheck,InterventionVector,param);
    InitialCond=Exsmooth(1,2);
    OriginalValues=100;
    MAxCondition=OriginalValues;
    [Mtrv]=RevSpaceTransform(Ncurve,InitialCond,100,25);
    DfferenceObs=MAxCondition(end)-Mtrv;
    % Update Initial State Values
    % Initial variance
    init_V(1,1)=max(R(2),param(2)^2);
    init_V(2,2)=param(4)^2*DfferenceObs+param(6)^2;
    init_V(3,3)=param(5)^2.;
end
% Filter
[Ex,VarKF,~,loglikely,MUSpeed,VarSpeed] = KalmanFilterFun(y, A, F, Q, Q_r, R, init_x, init_V,...
    ConstrainedKF,InterventionCheck,InterventionVector,...
    InterventionMu_Net,InterventionVar_Net,param);
% Smoother
TotalTimeSteps=length(Ex(1,:));
ExF(:,TotalTimeSteps)=Ex(:,TotalTimeSteps);
VarF(:,:,TotalTimeSteps)=VarKF(:,:,TotalTimeSteps);
[Exsmooth, Vsmooth, InterventionMu, InterventionVar,MUsSpeed,VarsSpeed]=KalmanSmootherFun(Ex,...
    VarKF,ExF,VarF,Q,Q_r,A,TotalTimeSteps,ConstrainedKF,InterventionCheck,...
    InterventionVector,param);
IntTime=find(InterventionVector);
if ~isempty(IntTime) && IntTime(end)>1
    for j=0:length(Exsmooth(1,IntTime(end):end))-1
        NCDF(j+1)=normcdf(Exsmooth(1,IntTime(end)-1),Exsmooth(1,IntTime(end)+j),sqrt(Vsmooth(1,1,IntTime(end)+j)));
    end
    NPDF=diff(NCDF);
    [~,IndvYear]=max(NPDF);
else
    NPDF=nan;
    NCDF=nan;
end
end

function [ExpectedCond, Variance, SigmaForecast, loglik,SpeedTracking_mu,SpeedTracking_var]=KalmanFilterFun(y,...
    A, F, Q, Q_r, R, init_x, init_V,ConstrainedKF,InterventionCheck,...
    InterventionVector,InterventionMu_Net,InterventionVar_Net,param)
T=length(y);
InType=1;
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
                ExpectedCond(4:6,t)=InterventionMu_Net{InType};
                Variance(4:6,4:6,t)=InterventionVar_Net{InType};
                Q_R=Q_r{InType};
                InType=InType+1;
            else
                A(1:3,4:6)=zeros(3);
                Q_R=zeros(6);
            end
        end
        UPvariance=A*Variance(:,:,t)*A'+ Q + Q_R;
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
        %%
        if InterventionVector(t+1)
            SpeedTracking_mu(1,1)=ExpectedCond(2,t+1);
            SpeedTracking_mu(2,1)=ExpectedCond(5,t+1);
            SpeedTracking_var(1,1)=Variance(2,2,t+1);
            SpeedTracking_var(2,1)=Variance(5,5,t+1);
        end
        %%
        KFConstraints();
        if InterventionVector(t+1)
            SpeedTracking_mu(3,1)=ExpectedCond(2,t+1);
            SpeedTracking_mu(4,1)=ExpectedCond(5,t+1);
            SpeedTracking_var(3,1)=Variance(2,2,t+1);
            SpeedTracking_var(4,1)=Variance(5,5,t+1);
        else
            SpeedTracking_mu=0;
            SpeedTracking_var=0;
        end
        SigmaForecast(:,t+1)=sqrt(diag(Variance(:,:,t+1)));
        if ~isnan(y(t+1))
            Err=y(t+1)-F*(A*ExpectedCond(:,t));
            loglik(t+1)=log(normpdf(Err,0,sqrt(S)));
        else
            loglik(t+1)=NaN;
        end
    end
end
loglik=sum(loglik(~isnan(loglik)));
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
InType=1;
for i=TotalTimeSteps-1:-1:Ending
    if length(InterventionVector)>1
        if InterventionVector(i+1)
            A(1:3,4:6)=eye(3);
            Q_R=Q_r{InType};
            InType=InType+1;
        else
            A(1:3,4:6)=zeros(3);
            Q_R=zeros(6);
        end
    end
    Xpred=A*x(:,i);
    Vpred=A*Var(:,:,i)*A' + Q + Q_R;
    J=Var(:,:,i)*A'*pinv(Vpred);
    ExSmooth(:,i)=x(:,i)+J*(ExSmooth(:,i+1)-Xpred);
    VarSmooth(:,:,i)=Var(:,:,i)+J*(VarSmooth(:,:,i+1)-Vpred)*J';
    ViIndex=find((2*sqrt(VarSmooth(2,2,i))+ExSmooth(2,i))>0);
    if InterventionVector(i+1)
        SpeedTracking_mu(1,1)=ExSmooth(2,i);
        SpeedTracking_mu(2,1)=ExSmooth(5,i);
        SpeedTracking_var(1,1)=VarSmooth(2,2,i);
        SpeedTracking_var(2,1)=VarSmooth(5,5,i);
    else
        SpeedTracking_mu=0;
        SpeedTracking_var=0;
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
    else
        InterventionMu=0;
        InterventionSigma=0;
    end

end
end