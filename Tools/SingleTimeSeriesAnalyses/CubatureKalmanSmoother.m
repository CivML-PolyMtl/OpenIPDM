%% Generic Code for Cubature Kalman Smoother
% Created by: Zachary Hamida
% November 19, 2019
%%
% This is a generic version of CKS, the original virsion is available
% at: https://www.haranarasaratnam.com/software.html
% Paper: Arasaratnam, Ienkaran, and Simon Haykin. "Cubature kalman smoothers."
% Automatica 47.10 (2011): 2245-2250.

%% INPUTS:
% 
% TransitionFunction:    non-linear transition function.
% Mu_KF: expected value of the state (Filter)
% Var_KF: variance of the state (Filter)
% Std_KF: standard deviation of the state (Filter)
% Q: process noise 

%% OUTPUTS:
%
% ExSmooth: expected value of the state (smoother)
% VarSmooth: variance of the state (smoother)

%% Note: This function is written to accommodate a single observation only
function [ExSmooth, VarSmooth]=CubatureKalmanSmoother(TransitionFunction, Mu_KF,...
    Var_KF,Std_KF,Q)
% size of state vector
nx=size(Mu_KF,1);
% Time series length
TotalTimeSteps=size(Mu_KF,2);
% Initilizing
Qsqrt=sqrt(Q);
[QPts,~,nPts]= FindCubaturePts(nx);
ExSmooth(:,TotalTimeSteps)=Mu_KF(:,TotalTimeSteps);
VarSmooth(:,:,TotalTimeSteps)=Var_KF(:,:,TotalTimeSteps);
StdSmooth(:,:,TotalTimeSteps)=Std_KF(:,:,TotalTimeSteps);

for t=TotalTimeSteps-1:-1:1
    X_sample = repmat(Mu_KF(:,t),1,nPts) + Std_KF(:,:,t)*QPts;
    Std_x_sample = (X_sample-repmat(Mu_KF(:,t),1,nPts))/sqrt(nPts);
    Var_x_sample = Std_x_sample*Std_x_sample';
    X_pred_sample = TransitionFunction(X_sample);
    Mu_Xpred_sample = sum(X_pred_sample,2)/nPts;
    Std_Xpred_sample = (X_pred_sample-repmat(Mu_Xpred_sample,1,nPts))/sqrt(nPts);
    Var_Xpred_sample = Std_Xpred_sample*Std_Xpred_sample';    
    Var_Xpred_update = Var_Xpred_sample + Q;
    J = Var_x_sample/Var_Xpred_update;
    ExSmooth(:,t) = Mu_KF(:,t) + J*( ExSmooth(:,t+1) - Mu_Xpred_sample );
    VarSmooth(:,:,t) = Var_KF(:,:,t) + J*(VarSmooth(:,:,t+1)-Var_Xpred_update)*J';
end
end

function [QPts,wPts,nPts]= FindCubaturePts(n)

nPts = 2*n;
wPts = ones(1, nPts)/(2*n);
QPts = sqrt(n)*eye(n);
QPts = [QPts -QPts];
end