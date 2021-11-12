%% Generic Code for Cubature Kalman Filter
% Created by: Zachary Hamida
% November 19, 2019
%%
% This is a generic version of CKF, the original virsion is available
% at: https://www.haranarasaratnam.com/software.html
% Paper: Arasaratnam, Ienkaran, and Simon Haykin. "Cubature kalman
% filters." IEEE Transactions on automatic control 54.6 (2009): 1254-1269.

%% INPUTS:
% 
% TransitionFunction:    Non-linear transition function.
% MeasurementFunction:   Non-linear Measurement function.
% y: Observation
% Q: Process noise
% R: Observation variance
% init_x: initial state at time t=0
% init_V: initial variance at time t=0

%% OUTPUTS:
%
% Mu_x: Expected value of the state
% Var_x: variance of the state
% Std_x: standard deviation of the state
% loglik: log likelihood of the observation

%% Note: This function is written to accommodate a single observation only


function [Mu_x, Var_x, Std_x, loglik]=CubtureKalmanFilter(TransitionFunction,...
    MeasurementFunction, y, Q, R, F, init_x, init_V)

% Calculate Cubature Points
nx=size(init_x,1);                       
[QPts,~,nPts]= FindCubaturePts(nx);

if ~isnan(y(1))
    y=[NaN y];
end

if length(R)<length(y)
    R=repmat(R,1,length(y));
end

% Time window
T=length(y)-1;

% Initilizatyion
loglik=zeros(T,1);

for t=0:T
    if t==0
        % Initilize the Kalman Filter
        Mu_x(:,1)=init_x;
        Var_x(:,:,1)=init_V;
        Std_x(:,:,1)=sqrt(Var_x);
    else
        % Sample from the current state
        X_sample = repmat(Mu_x(:,t),1,nPts) + Std_x(:,:,t)*QPts;
        % Propagate the sample in time
        X_sample = TransitionFunction(X_sample);
        % Mean and Std. of propgated sample
        Mu_x_sample = sum(X_sample,2)/nPts;
        Std_x_sample = (X_sample-repmat(Mu_x_sample,1,nPts))/sqrt(nPts);
        
        % Update the state variance with the process noise
        Var_x_update= Std_x_sample*Std_x_sample'+Q;
        [U, S, V] =  svd(Var_x_update);
        Std_x_update=0.5*(U+V)*sqrt(S);
        % Check for an observation
        if ~isnan(y(t+1))
            X_y_sample = MeasurementFunction(X_sample);
            Mu_x_y = sum(X_y_sample,2)/nPts;
            Std_x_y_sample = (X_y_sample-repmat(Mu_x_y,1,nPts))/sqrt(nPts);
            Var_x_y_sample = Std_x_y_sample*Std_x_y_sample';
            Var_x_y = Var_x_y_sample + R(t+1);
            Var_x_y_update=Var_x_update*F';
            % Compute Kalman gain
            G = Var_x_y_update/Var_x_y;
            % Update the mean and variance with the observation
            Mu_x(:,t+1) = Mu_x_sample + G*(y(t+1)-Mu_x_y);
            Var_x(:,:,t+1)=Var_x_update-G*Var_x_y*G';
            Std_x(:,:,t+1)=chol(Var_x(:,:,t+1),'lower');
        else
            Mu_x(:,t+1)=Mu_x_sample;
            Var_x(:,:,t+1)=Var_x_update;
            Std_x(:,:,t+1)=Std_x_update;
        end
    end
    if ~isnan(y(t+1))
        Err=y(t+1)-Mu_x_y;
        StdModel=sqrt(Var_x_y);
        loglik(t+1)=log(normpdf(Err,0,StdModel));
    else
        loglik(t+1)=NaN;
    end
end
loglik=sum(loglik(~isnan(loglik)));
end

function [QPts,wPts,nPts]= FindCubaturePts(n)
nPts = 2*n;
wPts = ones(1, nPts)/(2*n);
QPts = sqrt(n)*eye(n);
QPts = [QPts -QPts];
end