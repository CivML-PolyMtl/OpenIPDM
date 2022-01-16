function [x, V, VV, loglik,xinitc,Vinitc] = kalman_filter(y, A, C, Q, R, Re...
    , init_x, init_V, InpecBiase,BadInsp,RU,InspBU,ObsYears,Ncurve,OptBoundsData...
    ,ConstraintsKFs,SpeedConstraints,GpuCompute,Nsigma,varargin)
% Kalman filter.
% [x, V, VV, loglik] = kalman_filter(y, A, C, Q, R, init_x, init_V, ...)
%
% INPUTS:
% y(:,t)   - the observation at time t
% A - the system matrix
% C - the observation matrix 
% Q - the system covariance 
% R - the observation covariance
% init_x - the initial state (column) vector 
% init_V - the initial state covariance 
%
% OPTIONAL INPUTS (string/value pairs [default in brackets])
% 'model' - model(t)=m means use params from model m at time t [ones(1,T) ]
%     In this case, all the above matrices take an additional final dimension,
%     i.e., A(:,:,m), C(:,:,m), Q(:,:,m), R(:,:,m).
%     However, init_x and init_V are independent of model(1).
% 'u'     - u(:,t) the control signal at time t [ [] ]
% 'B'     - B(:,:,m) the input regression matrix for model m
%
% OUTPUTS (where X is the hidden state being estimated)
% x(:,t) = E[X(:,t) | y(:,1:t)]
% V(:,:,t) = Cov[X(:,t) | y(:,1:t)]
% VV(:,:,t) = Cov[X(:,t), X(:,t-1) | y(:,1:t)] t >= 2
% loglik = sum{t=1}^T log P(y(:,t))
%
% If an input signal is specified, we also condition on it:
% e.g., x(:,t) = E[X(:,t) | y(:,1:t), u(:, 1:t)]
% If a model sequence is specified, we also condition on it:
% e.g., x(:,t) = E[X(:,t) | y(:,1:t), u(:, 1:t), m(1:t)]


ss = size(A,1); % size of state space


u = [];
B = [];
ndx = [];

args = varargin;
nargs = length(args);
for i=1:2:nargs
  switch args{i}
   case 'model', model = args{i+1};
   case 'u', u = args{i+1};
   case 'B', B = args{i+1};
   case 'ndx', ndx = args{i+1};
   otherwise, error(['unrecognized argument ' args{i}])
  end
end

if GpuCompute
    [os Et T] = size(y);
    x = zeros(ss, Et, T,'gpuArray');
    V = zeros(ss, ss, Et, T,'gpuArray');
    VV = zeros(ss, ss, Et, T,'gpuArray');
else
    [os T] = size(y);
    x = zeros(ss, T);
    V = zeros(ss, ss, T);
    VV = zeros(ss, ss, T);
end
% set default params
model = ones(1,T);
loglik = 0;
for t=1:T

  m = model(t);
  if t==1
    %prevx = init_x(:,m);
    %prevV = init_V(:,:,m);
    prevx = init_x;
    prevV = init_V;
    initial = 1;
  else
      if GpuCompute
        prevx = x(:,:,t-1);
        prevV = V(:,:,:,t-1);
      else
          prevx = x(:,t-1);
          prevV = V(:,:,t-1);
      end
      initial = 0;
  end
  if isempty(u)
      if GpuCompute
        [x(:,:,t), V(:,:,:,t), LL, VV(:,:,:,t),xinitc(:,:,t),Vinitc(:,:,:,t)] = ...
        kalman_update(A(:,:,:,m), C(:,:,:,m), Q(:,:,:,m), R, Re(1,:,t), y(1,:,t), ...
        prevx, prevV, InpecBiase(1,:,t),BadInsp(:,t),RU(1,:,t),InspBU(1,:,t),ObsYears(1,:,t),...
        Ncurve,OptBoundsData,SpeedConstraints,ConstraintsKFs,GpuCompute,Nsigma,'initial', initial);
      else
          [x(:,t), V(:,:,t), LL, VV(:,:,t),xinitc(:,t),Vinitc(:,:,t)] = ...
        kalman_update(A(:,:,m), C(:,:,m), Q(:,:,m), R, Re(1,t), y(1,t), ...
        prevx, prevV, InpecBiase(1,t),BadInsp(t),RU(1,t),InspBU(1,t),ObsYears(1,t),...
        Ncurve,OptBoundsData,SpeedConstraints,ConstraintsKFs,GpuCompute,Nsigma,'initial', initial);
      end
  else
    if isempty(ndx)
      [x(:,t), V(:,:,t), LL, VV(:,:,t),xinitc(:,t),Vinitc(:,:,t)] = ...
	  kalman_update(A(:,:,m), C(:,:,m), Q(:,:,m), R(:,:,m), Re(t), y(:,t),...
      prevx, prevV, InpecBiase(t),BadInsp(t),RU(t),InspBU(t),ObsYears(t),Ncurve,...
      OptBoundsData,SpeedConstraints,ConstraintsKFs,Nsigma,'initial', initial, 'u', u(:,t), 'B', B(:,:,m));
    else
      i = ndx{t};
      % copy over all elements; only some will get updated
      x(:,t) = prevx;
      prevP = inv(prevV);
      prevPsmall = prevP(i,i);
      prevVsmall = inv(prevPsmall);
      [x(i,t), smallV, LL, VV(i,i,t),xinitc(:,t),Vinitc(:,:,t)] = ...
	  kalman_update(A(i,i,m), C(:,i,m), Q(i,i,m), R(:,:,m), Re(t), y(:,t), prevx(i),...
      prevVsmall, InpecBiase(t),BadInsp(t),RU(t),InspBU,ObsYears(t),Ncurve,OptBoundsData,...
      SpeedConstraints,ConstraintsKFs,Nsigma,...
			'initial', initial, 'u', u(:,t), 'B', B(i,:,m));
      smallP = inv(smallV);
      prevP(i,i) = smallP;
      V(:,:,t) = inv(prevP);
    end    
  end
  loglik = loglik + LL;
end



