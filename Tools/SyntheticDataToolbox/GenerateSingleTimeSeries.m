function [x,Y,RejectTimeSeries,RejectReason]=GenerateSingleTimeSeries(Condition,...
    Speed,Acc,FullTimeDuration,PARAM,Ncurve,MaxCond,MinCond,InspMaxStd,...
    InspMinStd,NumInsp,dt,InterventionsCheck,InterventionsParam,...
    MaxBias,MinBias)

SynInsp(1,:)=99001:(99000+NumInsp);                 % Inspector ID
SynInsp(2,:)=(InspMaxStd-InspMinStd).*rand(1,NumInsp) + InspMinStd;% Inspector Variance
if MinBias~=0 || MaxBias~=0
    SynInsp(3,:)=unifrnd(MinBias,MaxBias,1,NumInsp);                        % Inspector bias
else
    SynInsp(3,:)=0;
end
%% A - Model transition matrix
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];

%% Initial model parameters
Q=@(param) param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;(dt^4)/8 (dt^3)/3 (dt^2)/2;(dt^3)/6 (dt^2)/2 dt];
%% F - Observation Matrix
F=[1,0,0];      
% Generate engineers ID | engineers std and randomly assign groups of
% Repeat an inspector 4-5 times in each time-series
if rem(FullTimeDuration,2)>0
    InspRep=5;
else
    InspRep=4;
end
for i=1:1
    RepInsp=repelem(randperm(length(SynInsp(1,:)),FullTimeDuration/InspRep),1,InspRep);
    Shuffle1=SynInsp(:,RepInsp(randperm(FullTimeDuration)))'; %length(SynInsp(1,:))
    IdVarAll{i}=Shuffle1(1:FullTimeDuration,:);
end

% Interventions check
if InterventionsCheck
    % Interventions Parameters
    if isempty(InterventionsParam)
        load('InterventionsParam.mat');
    else
        In= InterventionsParam;
    end
    
    % load decision making system 
    DMS=readfis('DecisionMaker3.fis');
    InterventionTime=randi([5 9],1);
    %Priority: low: 0.5, mid:1.5, high: 2.5
    Priority=3.*rand(1,1);                
    Rehab_mu=[In.Rehab.Cond.mu;In.Rehab.Speed.mu;In.Rehab.Acc.mu];
    Rehab_cov=diag([In.Rehab.Cond.std;In.Rehab.Speed.std;In.Rehab.Acc.std]);
    PrM1_mu=[In.PrM1.Cond.mu;In.PrM1.Speed.mu;In.PrM1.Acc.mu];
    PrM1_cov=diag([In.PrM1.Cond.std;In.PrM1.Speed.std;In.PrM1.Acc.std]);
    PrM2_mu=[In.PrM2.Cond.mu;In.PrM2.Speed.mu;In.PrM2.Acc.mu];
    PrM2_cov=diag([In.PrM2.Cond.std;In.PrM2.Speed.std;In.PrM2.Acc.std]);
    % One time intervention
    InterventionDone=0;
end

% Generate Time-Series
[~,Ms]=SpaceTransformation(Ncurve,Condition,MaxCond,MinCond);
x0=[Ms;Speed;Acc];                             % random start state
LowerLimit=-4;
x(:,1)=x0;
x(:,60)=LowerLimit;
RejectTimeSeries=0;
RejectReason=0;
CountLoop=0;
CountLimit=200;

for t=1:FullTimeDuration
    InspectorVar=IdVarAll{1,1}(t,2);
    InspectorB=IdVarAll{1,1}(t,3);
    Constraint1=1;Constraint2=1;Constraint3=1;
    x(:,t+1)=(A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))');
    while x(2,t+1)>0&& x(3,t+1)>0.0001&& x(3,t+1)<-0.1...
            && Constraint1 && Constraint2 && Constraint3
        if x(2,t+1)>0
            x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
            CountLoop=CountLoop+1;
        else
            Constraint1=0;
            Constraint3=0;
            CountLoop=0;
        end
        if x(3,t+1)>0.0001
            x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
            CountLoop=CountLoop+1;
        else
            Constraint2=0;
            Constraint3=0;
            CountLoop=0;
        end
        if x(3,t+1)<-0.09
            x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
            CountLoop=CountLoop+1;
        else
            Constraint3=0;
            Constraint2=0;
            CountLoop=0;
        end
        if CountLoop>CountLimit
            Constraint1=0;
            Constraint2=0;
            Constraint3=0;
        end
    end
    if InterventionsCheck
        if t+1>=InterventionTime && ~InterventionDone
            RepCase=evalfis(DMS,[x(1,t+1), Priority]);

            if RepCase<=4 && RepCase>=3
                x(:,t+1)=x(:,t+1)+abs(mvnrnd(Rehab_mu,Rehab_cov)');
                InterventionDone=1;
            elseif RepCase<=2.999 && RepCase>=2
                x(:,t+1)=x(:,t+1)+abs(mvnrnd(PrM1_mu,PrM1_cov)');
                InterventionDone=1;
            else RepCase<=1.999 && RepCase>=1;
                x(:,t+1)=x(:,t+1)+abs(mvnrnd(PrM2_mu,PrM2_cov)');
                InterventionDone=1;
            end
        end
    end
    Y(t+1)=F*x(:,t+1)+normrnd(InspectorB,InspectorVar);
end
if isnan(x(1,end))
    RejectTimeSeries=1;
    RejectReason=1;
end
if CountLoop>CountLimit
    RejectTimeSeries=1;
    RejectReason=2;
end
if  x(2,1)<0.01*x(1,1)-1.3 %&& x(1,1)>=100*x(2,1)+110
    RejectTimeSeries=1;
    RejectReason=3;
end
if  ~isempty(find(x(2,:)>0))
    RejectTimeSeries=1;
    RejectReason=3;
end
% if  x(2,1)<-0.3 && x(1,1)>80
%     RejectTimeSeries=1;
%     RejectReason=4;
% end
% if  x(2,1)<-0.4 && x(1,1)>70
%     RejectTimeSeries=1;
%     RejectReason=5;
% end
if  x(3,1)<0.001*x(1,1)-0.13 %&& x(1,1)>=1000*x(3,1)+110
    RejectTimeSeries=1;
    RejectReason=6;
end
% if  x(3,1)<-0.03 && x(1,1)>80
%     RejectTimeSeries=1;
%     RejectReason=7;
% end
% if  x(3,1)<-0.04 && x(1,1)>70
%     RejectTimeSeries=1;
%     RejectReason=8;
% end
if x(1,1)<70
    RejectTimeSeries=1;
    RejectReason=9;
end
if x(1,end)>0.5*x(1,1) && ~InterventionsCheck
    RejectTimeSeries=1;
    RejectReason=10;
end
if x(1,FullTimeDuration/2)>0.9*x(1,1) && ~InterventionsCheck
    RejectTimeSeries=1;
    RejectReason=11;
end
if InterventionsCheck
    if x(1,end)>=x(1,InterventionTime)
        RejectTimeSeries=1;
        RejectReason=11;
    end
end
end