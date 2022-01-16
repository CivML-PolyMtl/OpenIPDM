function [ExpectedCond, Variance, SigmaForecast, loglik,SpeedTracking_mu,SpeedTracking_var]=KalmanFilterFun(y,...
    A, F, Q, Q_r, R, Bias, init_x, init_V,ConstrainedKF,InterventionCheck,...
    InterventionVector,InterventionMu_Net,InterventionVar_Net,param)
T=length(y);

% Equations
UpdateExpected=@(EX,EY,sigma_XY,sigma_obs,y,b) EX+sigma_XY*sigma_obs^-1*(y-EY-b);
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
                ExpectedCond(4:6,t)=InterventionMu_Net;
                Variance(4:6,4:6,t)=InterventionVar_Net;
                % added
%                 PredVal=A*ExpectedCond(:,t);
%                 OCond=RevSpaceTransform(4,PredVal(1));
%                 DCond=(100-OCond);
%                 Q_r([2,5],[2,5])=Q_r([2,5],[2,5]).*DCond+param^2*diag([1 1]);
                % added
                Q_R=Q_r;
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
            ExpectedCond(:,t+1)=UpdateExpected(A*ExpectedCond(:,t),F*(A*ExpectedCond(:,t)),sigmaXY,S,y(t+1),Bias(t+1)); 
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