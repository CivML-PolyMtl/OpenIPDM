% Input
SigmaWsliderValue = get(app.SigmaW,'Value');
InitialCondSigmasliderValue = get(app.InitialCondSigma,'Value');
InitialSpeedSigmasliderValue = get(app.InitialSpeedSigma,'Value');
InitialAccSigmasliderValue = get(app.InitialAccSigma,'Value');
ExpectedSpeedsliderValue = get(app.ExpectedSpeed,'Value');
ExpectedAccsliderValue = get(app.ExpectedAcc,'Value');
InterceptSpeedsliderValue = get(app.InterceptSpeed,'Value');
set(app.TxtSigmaW, 'Text',num2str(get(app.SigmaW,'Value')));
Ncurve=str2double(get(app.TransParam,'Value'));
% check intervention
InterventionCheck=get(app.UnderwentIntervention,'Value');
if InterventionCheck
    InterventionTime=get(app.InterventionSpinner,'Value');
    InterventionMu=get(app.IntmuEditField,'Value');
    InterventionSigma=get(app.IntsigmaEditField,'Value');
else
    InterventionTime=0;
    InterventionMu=0;
    InterventionSigma=0;
    InterventionVector=0;
end
set(app.TxtSigmaCond, 'Text',num2str(get(app.InitialCondSigma,'Value')));
set(app.TxtSigmaSpeed, 'Text',num2str(get(app.InitialSpeedSigma,'Value')));
set(app.TxtSigmaAcc, 'Text',num2str(get(app.InitialAccSigma,'Value')));
set(app.TxtExpectedSpeed, 'Text',num2str(get(app.ExpectedSpeed,'Value')));
set(app.TxtExpectedAcc, 'Text',num2str(get(app.ExpectedAcc,'Value')));
set(app.TxtIntercept, 'Text',num2str(get(app.InterceptSpeed,'Value')));

% Analysis Duration 
Duration=str2double(get(app.Durations,'Value'));

% observation set
DataTable=get(app.DataTable,'data');
y=[DataTable(:,1)' nan(1,Duration)];
R=[DataTable(:,2)' nan(1,Duration)].^2;

% transform values
% [~,y]=SpaceTransformation(Ncurve,y,100.01,25);

% Constraints
ConstrainedKF=1; % 1:Contrained, 0: No Constraints
if ~get(app.Constrained,'Value')
    ConstrainedKF=0;
end

% parameters case
Case1=get(app.Case1,'Value'); % default
Case2=get(app.Case2,'Value'); % other test
% parameters
% sigmaW
% sigma Init Cond
% sigma Init speed
% sigma Init Acc
% Expected Init Cond
% Expected Init speed
% Expected Init Acc
% Intercept (speed)
param=[SigmaWsliderValue
       InitialCondSigmasliderValue
       InitialSpeedSigmasliderValue
       InitialAccSigmasliderValue
       ExpectedSpeedsliderValue
       ExpectedAccsliderValue
       InterceptSpeedsliderValue
       ];

% Smoothed Plot vs. Non-Smoothed Plot
type=2; % 1:Filter only, 2: filter and smoother
originalspace=1; % 1: original space, 0: Transformed space

% time frame
timestamps=1:length(y);
timestamps=timestamps';

% Matrix A
dt=1;
Q=param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;(dt^4)/8 (dt^3)/3 (dt^2)/2;(dt^3)/6 (dt^2)/2 dt];
F=[1,0,0]; 

% Initial variance
if Case1
    A=[1 1 (dt^2)/2;0 1 dt;0 0 1];
    InitCase=1;
elseif Case2
    A=[1 0 0;0 1 0;0 0 1];
    Q=param(1).^2*[1.5983*(dt^5)/20 1.5983*(dt^4)/8 1.2642*(dt^3)/6;
                        1.5983*(dt^4)/8 1.5983*(dt^3)/3 1.2642*(dt^2)/2;
                        1.2642*(dt^3)/6 1.2642*(dt^2)/2 dt];
    InitCase=2;
end

% intervention check
if InterventionCheck
    A=blkdiag(A,eye(3));
    Q=blkdiag(Q,zeros(3));
    init_V=eye(6);
    init_V(4:6,4:6)=diag([InterventionSigma^2 10^-2 10^-2]);
    init_x=zeros(6,1);
    F=[1,zeros(1,5)];
    InterventionMu=[InterventionMu 0.1 0.01]';
    InterventionVector=zeros(1,length(y));
    InterventionVector(InterventionTime)=1;
end

% difference 
InitialCond=nanmean(y);%y(2);
OriginalValues=100;
MAxCondition=OriginalValues;
[Mtrv]=RevSpaceTransform(Ncurve,InitialCond);
% [~,MAxCondition]=SpaceTransformation(Ncurve,OriginalValues,100.01,25);
% InitialCond(find(InitialCond>MAxCondition(end)))=MAxCondition(end);
DfferenceObs=MAxCondition(end)-Mtrv;

% Initial expected value
init_x(1)=nanmean(y);%
init_x(2)=param(5)*DfferenceObs;
init_x(3)=param(6);%*init_x(2);
init_V(1,1)=max(R(2),param(2)^2);
init_V(2,2)=param(3)^2*DfferenceObs+param(7)^2;
init_V(3,3)=param(4)^2;


% KF - Kalman filter type
[Ex, Var, param, loglik, Exsmooth, Vsmooth]=HandCoded_KF(y,A,F,Q,R,param,...
    init_x,init_V,Ncurve,ConstrainedKF,InitCase,InterventionCheck,InterventionVector,...
    InterventionMu,InterventionSigma);   % Hand coded

% KF - Forecasts
ExF(:,1)=Ex(:,end);
VarF(:,:,1)=Var(:,:,end);

% for t=2:ForecastTime
%     ExF(:,t)=A*ExF(:,t-1);
%     VarF(:,:,t)=A*VarF(:,:,t-1)*A'+Q(param);
% end
for i=1:length(Exsmooth(1,:))
    Up2Std(1,i)=(RevSpaceTransform(Ncurve,2*sqrt(Vsmooth(1,1,i))+Exsmooth(1,i))-RevSpaceTransform(Ncurve,Exsmooth(1,i)));
    Down2Std(1,i)=(RevSpaceTransform(Ncurve,Exsmooth(1,i))-RevSpaceTransform(Ncurve,Exsmooth(1,i)-2*sqrt(Vsmooth(1,1,i))));
    ExKF(1,i)=(RevSpaceTransform(Ncurve,Exsmooth(1,i)));
    y(1,i)=RevSpaceTransform(Ncurve,y(1,i));
    Rtop(i)=RevSpaceTransform(Ncurve,2*sqrt(R(i))+y(i))-RevSpaceTransform(Ncurve,y(i));
    Rlow(i)=RevSpaceTransform(Ncurve,y(i))-RevSpaceTransform(Ncurve,y(i)-2*sqrt(R(i)));
end       
figure(1)
clf
figure(2)
clf
figure(3)
clf
if InterventionCheck
    PlotKF_Interventions
else
    PlotKF
end
