function [xnew, Vnew, loglik, VVnew,PlotSx2,PlotSy2,PlotSxy,PlotSxyg,xinitc,Vinitc] = kalman_update(A, C, Q, R, Re,y, x, V, InpecBiase,CurrentInspector,RU,InspBU,ObsYears,Ncurve,OptBoundsData,GlobalCondData,ConstraintKFs,Nsigma,varargin)
% KALMAN_UPDATE Do a one step update of the Kalman filter
% [xnew, Vnew, loglik] = kalman_update(A, C, Q, R, y, x, V, ...)
%
% INPUTS:
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% y(:)   - the observation at time t
% x(:) - E[X | y(:, 1:t-1)] prior mean
% V(:,:) - Cov[X | y(:, 1:t-1)] prior covariance
%
% OPTIONAL INPUTS (string/value pairs [default in brackets])
% 'initial' - 1 means x and V are taken as initial conditions (so A and Q are ignored) [0]
% 'u'     - u(:) the control signal at time t [ [] ]
% 'B'     - the input regression matrix
%
% OUTPUTS (where X is the hidden state being estimated)
%  xnew(:) =   E[ X | y(:, 1:t) ] 
%  Vnew(:,:) = Var[ X(t) | y(:, 1:t) ]
%  VVnew(:,:) = Cov[ X(t), X(t-1) | y(:, 1:t) ]
%  loglik = log P(y(:,t) | y(:,1:t-1)) log-likelihood of innovatio

% set default params
u = [];
B = [];
initial = 0;
xinitc=zeros(3,1);
Vinitc=zeros(3,3);
args = varargin;
for i=1:2:length(args)
  switch args{i}
   case 'u', u = args{i+1};
   case 'B', B = args{i+1};
   case 'initial', initial = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

%  xpred(:) = E[X_t+1 | y(:, 1:t)]
%  Vpred(:,:) = Cov[X_t+1 | y(:, 1:t)]

if initial
    
  if isempty(u)
    xpred = x;
  else
    xpred = x + B*u;
  end
  Vpred = V;
  
else
  if isempty(u)
    xpred = A*x;
  else
    xpred = A*x + B*u;
  end
  Vpred = A*V*A' + Q;
end
if ObsYears==0
    xnew=xpred;
    Vnew=Vpred;
    loglik=0;
    VVnew=A*V;
    PlotSx2=Vpred;
    S = C*Vpred*C' + R;
    PlotSy2=sqrt(diag(Vpred));%S;
    PlotSxy=Vpred*C';
    PlotSxyg=Vpred-PlotSxy*S^-1*PlotSxy';
else
    if CurrentInspector==1
        e = y - C*xpred - InspBU; % error (innovation) 
        R=RU;
    elseif CurrentInspector==0
        e = y - C*xpred - InpecBiase; % error (innovation)
        R=Re;
    else
        e = y - C*xpred; % error (innovation)
        R=Re;
    end
    n = length(e);
    ss = length(A);
    S = C*Vpred*C' + R;
    Sinv = inv(S);
    ss = length(V);
    loglik = gaussian_prob(e, zeros(1,length(e)), S, 1);
    K = Vpred*C'*Sinv; % Kalman gain matrix
    % If there is no observation vector, set K = zeros(ss).
    xnew = xpred + K*e;
    Vnew = (eye(ss) - K*C)*Vpred;
    VVnew = (eye(ss) - K*C)*A*V;
    PlotSx2=Vpred;
    PlotSy2=sqrt(diag(Vpred));%S;
    PlotSxy=Vpred*C';
    PlotSxyg=Vpred-PlotSxy*S^-1*PlotSxy';
end
    % Added section
%     if xnew(1)>100 || xnew(1)<25 || xnew(2)>0 || xnew(2)<-50
    if ConstraintKFs && ~initial
        if (2*sqrt(Vnew(2,2))+xnew(2))>GlobalCondData(2,4)
            d=[GlobalCondData(2,3);GlobalCondData(2,4)];
            D=[0 1 0;0 1 0];
            [xnew,Vnew,Status]=KFConstraintsHandling(xnew,Vnew,D,d,Nsigma);
        end
        if (-2*sqrt(Vnew(2,2))+xnew(2))<GlobalCondData(2,3)
            d=[GlobalCondData(2,3);GlobalCondData(2,4)];
            D=[0 1 0;0 1 0];
            [xnew,Vnew,Status]=KFConstraintsHandling(xnew,Vnew,D,d,Nsigma);
        end
    end
    if isnan(loglik)
        loglik=-100000;
    end
end
