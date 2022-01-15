function [xnew, Vnew, loglik, VVnew,xinitc,Vinitc] = kalman_update(A, C, Q...
    , R, Re,y, x, V, InpecBiase,CurrentInspector,RU,InspBU,ObsYears,Ncurve...
    ,OptBoundsData,GlobalCondData,ConstraintKFs,GPUCompute,Nsigma,varargin)
% KALMAN_UPDATE Do a one step update of the Kalman filter
% [xnew, Vnew, loglik] = kalman_update(A, C, Q, R, y, x, V, ...)


% set default params
u = [];
B = [];
initial = 0;

if GPUCompute
    xinitc=zeros(3,length(x));
    Vinitc=zeros(3,3,length(x));
else
    xinitc=zeros(length(x),1);
    Vinitc=zeros(length(x));
end
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
      if GPUCompute
        xpred = pagefun(@mtimes,A,x);
      else
        xpred = A*x;
      end
  else
    xpred = A*x + B*u;
  end
  if GPUCompute
    Vpred = pagefun(@mtimes,pagefun(@mtimes,A,V),pagefun(@transpose,A)) + Q;
  else
    Vpred = A*V*A' + Q;  
  end
end
if GPUCompute
    loglik=zeros(length(ObsYears(1,:,1)),1,'gpuArray');
    e=NaN(length(ObsYears(1,:,1)),1,'gpuArray');
    xnew=xpred;
    Vnew=Vpred;
    VVnew=pagefun(@mtimes,A,V);
else
    loglik=0;
    e=NaN(length(ObsYears),1);
    xnew=xpred;
    Vnew=Vpred;
    VVnew=A*V;
end

InspInd2=find(~isnan(y));   
if ~isempty(InspInd2)
    if GPUCompute
        e(InspInd2) = y(InspInd2) - pagefun(@mtimes,C,xpred(:,InspInd2))-InpecBiase(InspInd2); % error (innovation) 
        R(InspInd2)= Re(InspInd2);
    else
        e = y - C*xpred-InpecBiase; % error (innovation) 
        R= Re;
    end
end

InspInd=find(CurrentInspector==1);
if ~isempty(InspInd)
    if GPUCompute
        e(InspInd) = y(InspInd) - pagefun(@mtimes,C,xpred(:,InspInd))-InspBU(InspInd); % error (innovation) 
        R(InspInd)= RU(InspInd);
    else
        e = y - C*xpred-InspBU; % error (innovation) 
        R= RU;
    end
end

UpdateInd=union(InspInd,InspInd2);
if ~isempty(UpdateInd)
    n = length(e(UpdateInd));
    ss = size(A,1);
    if GPUCompute
        S = pagefun(@mtimes,pagefun(@mtimes,C,Vpred(:,:,UpdateInd)),pagefun(@transpose,C)) + reshape(R(UpdateInd),[],1,length(R(UpdateInd)));
        Sinv = pagefun(@inv,S);
        LikeTr=1;
        loglik(UpdateInd)=log(normpdf(e(UpdateInd), 0, sqrt(reshape(S,length(S),1,[])))).*LikeTr';
         K = pagefun(@mtimes,pagefun(@mtimes,Vpred(:,:,UpdateInd),pagefun(@transpose,C)),Sinv); % Kalman gain matrix
        % If there is no observation vector, set K = zeros(ss).
        xnew(:,UpdateInd) = xpred(:,UpdateInd) + reshape(pagefun(@mtimes,K,reshape(e(UpdateInd),[],1,length(e(UpdateInd)))),3,length(UpdateInd));
        Vnew(:,:,UpdateInd) = pagefun(@mtimes,(eye(ss) - pagefun(@mtimes,K,C)),Vpred(:,:,UpdateInd));
        VVnew(:,:,UpdateInd) = pagefun(@mtimes,(eye(ss) - pagefun(@mtimes,K,C)),pagefun(@mtimes,A,V(:,:,UpdateInd)));
    else
        S = C*Vpred*C' + R;
        Sinv =inv(S);
        loglik=log(normpdf(e, 0, sqrt(S)));
        K = Vpred*C'*Sinv; % Kalman gain matrix
        xnew = xpred + K*e;
        Vnew = (eye(ss) - K*C)*Vpred;
        VVnew = (eye(ss) - K*C)*A*V;
    end
end
    % Added section
%     if xnew(1)>100 || xnew(1)<25 || xnew(2)>0 || xnew(2)<-50
    if ConstraintKFs && ~initial
        ViIndex=find((2*sqrt(Vnew(2,2,:))+reshape(xnew(2,:),1,1,length(xnew(2,:))))>GlobalCondData(2,4));
    if ~isempty(ViIndex)
        xnewIN=xnew(:,ViIndex);
        VnewIN=Vnew(:,:,ViIndex);
        d=[GlobalCondData(2,3);GlobalCondData(2,4)];
        D=[0 1 0;0 1 0];
        if GPUCompute
            [xnewIN,VnewIN,Status]=KFConstraintsHandling(xnewIN,VnewIN,D,d,Nsigma);
            xnew(:,ViIndex)=xnewIN;
            Vnew(:,:,ViIndex)=VnewIN;
        else
            [xnewIN,VnewIN,Status]=KFConstraintsHandlingSeq(xnewIN,VnewIN,D,d,Nsigma);
            xnew=xnewIN;
            Vnew=VnewIN;
        end
    end
        ViIndex=find((-2*sqrt(Vnew(2,2,:))+reshape(xnew(2,:),1,1,length(xnew(2,:))))<GlobalCondData(2,3));
    if ~isempty(ViIndex)
        xnewIN=xnew(:,ViIndex);
        VnewIN=Vnew(:,:,ViIndex);
        d=[GlobalCondData(2,3);GlobalCondData(2,4)];
        D=[0 1 0;0 1 0];
        if GPUCompute
            [xnewIN,VnewIN,Status]=KFConstraintsHandling(xnewIN,VnewIN,D,d,Nsigma);
            xnew(:,ViIndex)=xnewIN;
            Vnew(:,:,ViIndex)=VnewIN;
        else
            [xnewIN,VnewIN,Status]=KFConstraintsHandlingSeq(xnewIN,VnewIN,D,d,Nsigma);
            xnew=xnewIN;
            Vnew=VnewIN;
        end
    end
    end
    NanlogIndex=find(isnan(loglik(UpdateInd)));
    if ~isempty(NanlogIndex)
        loglik(UpdateInd(NanlogIndex))=-1000000;
    end
end
