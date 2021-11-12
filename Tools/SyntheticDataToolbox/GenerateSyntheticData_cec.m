function [TrueState,ObservedData,GeneratedSyntheticData]=...
    GenerateSyntheticData_cec(app,T,dt,DatasetSize,NumInsp,InspMinStd,...
    InspMaxStd,PARAM,Sigma_w1,Sigma_w2,MaxStart,MaxTimeCut,StdInitSpeedV,EInitSpeedV,...
    StdInitAccV,EInitAccV,Ncurve,MaxCond,MinCond,InterventionsCheck,WeightedProb)

% Number of Observations per time series: weights
if isempty(WeightedProb)
    WeightedProb=[0,0,0.5403,0.3412,0.0940,0.0244,1.1998e-04,0,0,1.9980e-05];
end
%% Time Series Length
TimeSeriesWindow=length(WeightedProb);

%% Inpsector Data
rng(1,'twister');
SynInsp(1,:)=99001:(99000+NumInsp);                                         % Inspector ID
SynInsp(2,:)=(InspMaxStd-InspMinStd).*rand(1,NumInsp) + InspMinStd;         % Inspector Variance


%% A - transition matrix
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];

%% Full duration
FullTimeDuration=T+MaxTimeCut;
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
for i=1:DatasetSize
    RepInsp=repelem(randperm(length(SynInsp(1,:)),round(FullTimeDuration/InspRep)),1,InspRep);
    Shuffle1=SynInsp(:,RepInsp(randperm(FullTimeDuration)))'; %length(SynInsp(1,:))
    IdVarAll{i}=Shuffle1(1:FullTimeDuration,:);
end

% sample a starting point
% propotions
if InterventionsCheck
    % LowestCond: represent the lowest condition to start time series
    Cond_Percent=[0.1 0.3 0.6]; LowestCond=50; MidCond=70;
else
    Cond_Percent=[0.3 0.6 0.1]; LowestCond=70; MidCond=80;
end
% counters to ensure diversity
Initial90_100=round(DatasetSize*Cond_Percent(1));
Initial80_90=round(DatasetSize*Cond_Percent(2));
Initial70_80=round(DatasetSize*Cond_Percent(3));
Count90_100=0;Count80_90=0;Count70_80=0;

MinTimeCut=40;
SplitWindow=round((MaxTimeCut-MinTimeCut)/3);


% Interventions check
if InterventionsCheck
    % Interventions Parameters
    if isempty(app.InterventionsParam)
        load('InterventionsParam.mat');
    else
        In=app.InterventionsParam;
    end
    
    % load decision making system 
    DMS=readfis('DecisionMaker3.fis');
    InterventionTime_0=randi([5 9],DatasetSize,1);
    %Priority: low: 0, mid:1.5, high: 3
    Priority=3.*rand(1,DatasetSize);             
    Rehab_mu=[In.Rehab.Cond.mu;In.Rehab.Speed.mu;In.Rehab.Acc.mu];
    Rehab_cov=diag([In.Rehab.Cond.std;In.Rehab.Speed.std;In.Rehab.Acc.std].^2);
    PrM1_mu=[In.PrM1.Cond.mu;In.PrM1.Speed.mu;In.PrM1.Acc.mu];
    PrM1_cov=diag([In.PrM1.Cond.std;In.PrM1.Speed.std;In.PrM1.Acc.std].^2);
    PrM2_mu=[In.PrM2.Cond.mu;In.PrM2.Speed.mu;In.PrM2.Acc.mu];
    PrM2_cov=diag([In.PrM2.Cond.std;In.PrM2.Speed.std;In.PrM2.Acc.std].^2);
    
    % one time intervention
    InterventionDone=0;
end

SynDatabase={};
InspecData=[];InspecDataShort=[];
[~,MaxVCondition]=SpaceTransformation(Ncurve,100,MaxCond,MinCond);
CountLoop=0;
CountLimit=200;

for i=1:DatasetSize
    % speed initial state
    SsSigma_ln=sqrt(log(1+(StdInitSpeedV/abs(EInitSpeedV))^2));
    SsExpected_ln=log(abs(EInitSpeedV))-0.5*SsSigma_ln^2;
    % sample from initial speed state
    Ss=-lognrnd(SsExpected_ln,SsSigma_ln);
    % acceleration initial state
    AsSigma_ln=sqrt(log(1+(StdInitAccV/abs(EInitAccV))^2));
    AsExpected_ln=log(abs(EInitAccV))-0.5*AsSigma_ln^2;
    % sample from initial acceleration state
    As=-lognrnd(AsExpected_ln,AsSigma_ln);
    % transformed space condition
    [~,Ms]=SpaceTransformation(Ncurve,MaxStart,MaxCond,MinCond);
    x0=[Ms;Ss;As];
    LowerLimit=-4;                                                         
    x(:,1)=x0;
    x(:,60)=LowerLimit;
    RejectTimeSeries=1;
    InterventionDone=0;
    clc
    (fprintf('Number of accepted time series %d/%d',i,DatasetSize));
    while RejectTimeSeries
        RejectTimeSeries=0;
        WeightedProbInd=1-[repmat(Count90_100/Initial90_100,1,SplitWindow)...
            repmat(Count80_90/Initial80_90,1,SplitWindow) ...
            repmat(Count70_80/Initial70_80,1,MaxTimeCut-2*SplitWindow)];
        TimeCut=randsample(MinTimeCut:MaxTimeCut,1,true,WeightedProbInd);
        % Observations sample
        SampledIDs=sort(randperm(TimeSeriesWindow,randsample(TimeSeriesWindow,1,true,WeightedProb)));
        for t=1:FullTimeDuration
            InspectorVar=IdVarAll{1,i}(t,2);
            Constraint1=1;Constraint2=1;Constraint3=1;
            x(:,t+1)=(A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))');
            
            
            while x(2,t+1)>0&& x(3,t+1)>0.0001&& x(3,t+1)<-0.1 && Constraint1 && Constraint2 && Constraint3
                % check positive speed
                if x(2,t+1)>0
                    x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
                    CountLoop=CountLoop+1;
                else
                    Constraint1=0;
                    Constraint3=0;
                    CountLoop=0;
                end
                % check positive acceleration
                if x(3,t+1)>0.0001
                    x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
                    CountLoop=CountLoop+1;
                else
                    Constraint2=0;
                    Constraint3=0;
                    CountLoop=0;
                end
                % check highly negative acceleration
                if x(3,t+1)<-0.09
                    x(:,t+1)=A*x(:,t)+mvnrnd(zeros(1,3),Q(PARAM))';
                    CountLoop=CountLoop+1;
                else
                    Constraint3=0;
                    Constraint2=0;
                    CountLoop=0;
                end
                % check count limit (ensure diversity)
                if CountLoop>CountLimit
                    Constraint1=0;
                    Constraint2=0;
                    Constraint3=0;
                end
                
            end
            if InterventionsCheck
                InterventionTime=InterventionTime_0(i)+TimeCut;
                if t+1==InterventionTime && ~InterventionDone
                    CondOriginal=RevSpaceTransform(Ncurve,x(1,t+1),MaxCond,MinCond);
                    RepCase=evalfis(DMS,[CondOriginal, Priority(i)]);
                    % intervention catigoery
                    if (RepCase<=4 && RepCase>=3)
                        x(:,t+1) = x(:,t+1)+ Rehab_mu + mvnrnd(zeros(3,1),Rehab_cov)';
                        InterventionDone=1;
                        InterventionCase(i)=3;
                        
                    elseif (RepCase<=2.999 && RepCase>=2)
                        x(:,t+1) = x(:,t+1)+ PrM1_mu + mvnrnd(zeros(3,1),PrM1_cov)';
                        InterventionDone=1;
                        InterventionCase(i)=2;
                    elseif RepCase<=1.999 && RepCase>=1
                        x(:,t+1) = x(:,t+1)+ PrM2_mu + mvnrnd(zeros(3,1),PrM2_cov)';
                        InterventionDone=1;
                        InterventionCase(i)=1;
                    end
%                     InterventionState(:,i)=x(:,t);
%                     %InterventionTime_0(i)=t+1-TimeCut;
                end
            end
            Y(t+1)=F*x(:,t+1)+normrnd(0,InspectorVar);
            Y(t+1)=RevSpaceTransform(Ncurve,Y(t+1),MaxCond,MinCond);
            if Y(t+1)>=90
                Y(t+1)=95;
            elseif Y(t+1)<90 && Y(t+1)>=80
                Y(t+1)=85;
            elseif Y(t+1)<80 && Y(t+1)>=70   
                Y(t+1)=75;
            else
                Y(t+1)=35;
            end
        end
        
%         TimeCut=randi([30,MaxTimeCut],1,1);
        % trim the time-series
        x(:,1:TimeCut)=[];
        x(:,T+1:end)=[];
        % check nan
        if isnan(x(1,end))
            RejectTimeSeries=1;
        end
        % check count (ensure diversity)
        if CountLoop>CountLimit
            RejectTimeSeries=1;
        end
        % speed and acceleration check
        if  x(2,1)<0.01*x(1,1)-1.3 %&& x(1,1)>=100*x(2,1)+110
            RejectTimeSeries=1;
        end
        if  x(3,1)<0.001*x(1,1)-0.13 %&& x(1,1)>=1000*x(3,1)+110
            RejectTimeSeries=1;
        end
        % check positive speed 
        if  x(2,1)>0
            RejectTimeSeries=1;
        end
        if  ~isempty(find(x(2,:)>0))
            RejectTimeSeries=1;
        end
        
%         if  x(2,1)<0.01*x(1,1)-1.1 && x(1,1)>=100*x(2,1)-110
%             RejectTimeSeries=1;
%         end
%         if  x(2,1)<-0.3 && x(1,1)>80
%             RejectTimeSeries=1;
%         end
%         if  x(2,1)<-0.4 && x(1,1)>70
%             RejectTimeSeries=1;
%         end
%         if  x(3,1)<-0.02 && x(1,1)>90
%             RejectTimeSeries=1;
%         end
%         if  x(3,1)<-0.03 && x(1,1)>80
%             RejectTimeSeries=1;
%         end
%         if  x(3,1)<-0.04 && x(1,1)>70
%             RejectTimeSeries=1;
%         end
        % check lowest possible condition at a time step
        if x(1,1)<LowestCond % 70
            RejectTimeSeries=1;
        end
        if x(1,end)>0.5*x(1,1) && ~InterventionsCheck
            RejectTimeSeries=1;
        end
        if InterventionDone
            if x(1,end)>0.75*x(1,InterventionTime_0(i)) && InterventionsCheck 
                RejectTimeSeries=1;
            end
        end
        if x(1,round(T/2))>0.85*x(1,1) && ~InterventionsCheck
            RejectTimeSeries=1;
        end
        if InterventionDone
            if x(1,round(T/1.5))>0.85*x(1,InterventionTime_0(i)) && InterventionsCheck 
                RejectTimeSeries=1;
            end
        end
        
        % check interventions
        if InterventionsCheck
            if ~InterventionDone
                RejectTimeSeries=1;
            else
                if x(1,end)>=x(1,InterventionTime_0(i))
                    RejectTimeSeries=1;
                end
                if x(2,InterventionTime_0(i)-1)>x(2,InterventionTime_0(i))
                    RejectTimeSeries=1;
                end
                if x(3,InterventionTime_0(i)-1)>=0
                    RejectTimeSeries=1;
                end
                if sum((SampledIDs>InterventionTime_0(i)))>=2 && length(SampledIDs)<=3
                    RejectTimeSeries=1;
                elseif sum((SampledIDs>InterventionTime_0(i)))>1 && length(SampledIDs)<3
                    RejectTimeSeries=1;
                elseif sum((SampledIDs>InterventionTime_0(i)))==0
                    RejectTimeSeries=1;
                end
            end
        end
        
        if x(1,1)>90 && RejectTimeSeries==0
            if Count90_100>=Initial90_100
                RejectTimeSeries=1;
            else
                Count90_100=Count90_100+1;
            end
        elseif x(1,1)<=90 && x(1,1)>MidCond && RejectTimeSeries==0
            if Count80_90>=Initial80_90
                RejectTimeSeries=1;
            else
                Count80_90=Count80_90+1;
            end
        elseif x(1,1)<=MidCond && x(1,1)>=LowestCond && RejectTimeSeries==0
            if Count70_80>=Initial70_80
                RejectTimeSeries=1;
            else
                Count70_80=Count70_80+1;
            end
        end
        if RejectTimeSeries
            Ss=-lognrnd(SsExpected_ln,SsSigma_ln);
%             Ss=unifrnd(EInitSpeedV,StdInitSpeedV);
            As=-lognrnd(AsExpected_ln,AsSigma_ln);
%             As=unifrnd(EInitAccV,StdInitAccV);
%             As=normrnd(EInitAccV,StdInitAccV);
            x0=[Ms;Ss;As];
            x(:,1)=x0;
            x(1,end)=LowerLimit;
            CountLoop=0;
            if InterventionsCheck
                % one time intervention
                InterventionDone=0;
            end
        end
    end
    Y(1:TimeCut+1)=[];
    Y(T:end)=[];
    IdVarAll{1,i}(1:TimeCut,:)=[];
    IdVarAll{1,i}(T:end,:)=[];
    if app.GenerateStructuralAttributesDefault2AttributesCheckBox.Value && InterventionsCheck
        StracturalAttributeA(i)=log(abs(x(2,SampledIDs(1))))+normrnd(0,Sigma_w1);
        StracturalAttributeB(i)=exp(x(2,SampledIDs(1)))+normrnd(0,Sigma_w2);
        SynDatabase{i,1}=[Y' IdVarAll{1,i}(:,1) [1965:(1964+T-1)]' i*ones(T-1,1) StracturalAttributeB(i)*ones(T-1,1) StracturalAttributeA(i)*ones(T-1,1) (InterventionTime_0(i)+1964-1)*ones(T-1,1) InterventionCase(i)*ones(T-1,1)];
    elseif app.GenerateStructuralAttributesDefault2AttributesCheckBox.Value && ~InterventionsCheck
        StracturalAttributeA(i)=log(abs(x(2,SampledIDs(1))))+normrnd(0,Sigma_w1);
        StracturalAttributeB(i)=exp(x(2,SampledIDs(1)))+normrnd(0,Sigma_w2);
        SynDatabase{i,1}=[Y' IdVarAll{1,i}(:,1) [1965:(1964+T-1)]' i*ones(T-1,1) StracturalAttributeB(i)*ones(T-1,1) StracturalAttributeA(i)*ones(T-1,1)];
    elseif ~app.GenerateStructuralAttributesDefault2AttributesCheckBox.Value && InterventionsCheck
        SynDatabase{i,1}=[Y' IdVarAll{1,i}(:,1) [1965:(1964+T-1)]' i*ones(T-1,1) (InterventionTime_0(i)+1964-1)*ones(T-1,1) InterventionCase(i)*ones(T-1,1)];
    else
        SynDatabase{i,1}=[Y' IdVarAll{1,i}(:,1) [1965:(1964+T-1)]' i*ones(T-1,1)];
    end
    
%     SampledIDs=sort(randperm(10,randi([3,6],1,1)));
%     WeightedProb=[0,0,0.77,0.1815,0.0472,6.03E-4,0,6.97E-4,0,0];
     
    SynDatabaseShort{i,1}=SynDatabase{i,1}(SampledIDs,:);
    InspecData=[InspecData IdVarAll{1,i}(:,1)];
    InspecDataShort{i,1}=IdVarAll{1,i}(SampledIDs,1);
    InitialX2(i)=x(2,1);
    InitialX2v(i)=x(2,SampledIDs(1));
    SynDatabaseState{i}=x;
    if InterventionsCheck
        InterventionState(:,i)=x(:,InterventionTime_0(i)-1);
        InterventionCategoryRec(i) = InterventionCase(i);
    end
%     SumY=[Y' SumY];
    x=[];
    Y=[];
    ySpeedTransformed(i)=SynDatabaseState{1,i}(2,SampledIDs(1));
    yAccTransformed(i)=SynDatabaseState{1,i}(3,SampledIDs(1));
    xObsTransformedAv(i)=MaxVCondition-SynDatabaseShort{i,1}(1,1);%SynDatabaseState{1,i}(1,SampledIDs(1));%
    xObsTransformedTrue(i)=MaxVCondition-SynDatabaseState{1,i}(1,SampledIDs(1)+1);%
    xSpeedTransformedTrue(i)=SynDatabaseState{1,i}(2,SampledIDs(1));%
end
figure(3)
plot(xObsTransformedAv,ySpeedTransformed,'o');
hold on
coeffs = polyfit(xObsTransformedAv, ySpeedTransformed, 1);
% Get fitted values
fittedX = linspace(min(xObsTransformedAv), max(xObsTransformedAv), 200);
fittedY = polyval(coeffs, fittedX);
% Plot the fitted line
plot(fittedX, fittedY, 'LineWidth', 2);

figure(3)
plot(xObsTransformedTrue,ySpeedTransformed,'o');
hold on
coeffs1 = polyfit(xObsTransformedTrue, ySpeedTransformed, 1);
% Get fitted values
fittedX = linspace(min(xObsTransformedTrue), max(xObsTransformedTrue), 200);
fittedY = polyval(coeffs1, fittedX);
% Plot the fitted line
plot(fittedX, fittedY,'-b', 'LineWidth', 3);
legend('Observation','Line Fit','True Condition','True Line Fit');

figure(4)
plot(xObsTransformedAv,yAccTransformed,'o');
hold on
coeffa = polyfit(xObsTransformedAv, yAccTransformed, 1);
% Get fitted values
fittedX = linspace(min(xObsTransformedAv), max(xObsTransformedAv), 200);
fittedY = polyval(coeffa, fittedX);
% Plot the fitted line
plot(fittedX, fittedY, 'LineWidth', 2);

figure(4)
plot(xObsTransformedTrue,yAccTransformed,'o');
hold on
coeffa1 = polyfit(xObsTransformedTrue, yAccTransformed, 1);
% Get fitted values
fittedX = linspace(min(xObsTransformedTrue), max(xObsTransformedTrue), 200);
fittedY = polyval(coeffa1, fittedX);
% Plot the fitted line
plot(fittedX, fittedY,'-b', 'LineWidth', 3);
legend('Observation','Line Fit','True Condition','True Line Fit');

figure(5)
plot(xSpeedTransformedTrue,yAccTransformed,'o');
hold on
coeffa = polyfit(xSpeedTransformedTrue, yAccTransformed, 1);
% Get fitted values
fittedX = linspace(min(xSpeedTransformedTrue), max(xSpeedTransformedTrue), 200);
fittedY = polyval(coeffa, fittedX);
% Plot the fitted line
plot(fittedX, fittedY, 'LineWidth', 2);
legend('True Speed vs. True Acc.','True Line Fit');

if app.GenerateStructuralAttributesDefault2AttributesCheckBox.Value
    figure(7)
    plot(InitialX2,StracturalAttributeA,'o');
    xlabel('Initial Speed')
    ylabel('Structural Attribute A')
    figure(8)
    plot(InitialX2,StracturalAttributeB,'o');
    xlabel('Initial Speed')
    ylabel('Structural Attribute B')
    
end

SampleID=randi([1,DatasetSize],1,100);
TrueState=SynDatabaseState(1,SampleID);
ObservedData=SynDatabase(SampleID,1);
% TrueState=SynDatabaseState(1,:);
% ObservedData=SynDatabase(:,1);
GeneratedSyntheticData.TrueState=SynDatabaseState';
GeneratedSyntheticData.FullObserved=SynDatabase;
GeneratedSyntheticData.ShortObserved=SynDatabaseShort;
GeneratedSyntheticData.InspectorsData=InspecData;
GeneratedSyntheticData.InspectorsDataShort=InspecDataShort;
GeneratedSyntheticData.InspectorsDataError=SynInsp';
GeneratedSyntheticData.SigmaW=PARAM;
GeneratedSyntheticData.TimeSpan=T;
GeneratedSyntheticData.TimeCut=MaxTimeCut;
GeneratedSyntheticData.TimeStep=dt;
GeneratedSyntheticData.Speed_Init=[EInitSpeedV,StdInitSpeedV];
GeneratedSyntheticData.Speed_CoefObs=coeffs;
GeneratedSyntheticData.Speed_CoefTrue=coeffs1;
GeneratedSyntheticData.Acc_Init=[EInitAccV,StdInitAccV];
GeneratedSyntheticData.Acc_CoefObs=coeffa;
GeneratedSyntheticData.Acc_CoefTrue=coeffa1;
GeneratedSyntheticData.NCurveparam=Ncurve;
if app.GenerateStructuralAttributesDefault2AttributesCheckBox.Value
    GeneratedSyntheticData.StracturalAttributeA=StracturalAttributeA;
    GeneratedSyntheticData.StracturalAttributeB=StracturalAttributeB;
    GeneratedSyntheticData.StracturalAttributeANoise=Sigma_w1;
    GeneratedSyntheticData.StracturalAttributeBNoise=Sigma_w2;
end
if InterventionsCheck
    GeneratedSyntheticData.InterventionCase=InterventionCase;
    GeneratedSyntheticData.InterventionTime=InterventionTime_0+1964;
    GeneratedSyntheticData.InterventionParam=In;
    GeneratedSyntheticData.PriorityMetric=Priority;
    GeneratedSyntheticData.DecisionMaker=DMS;
end
% set(handles.SyntheticDataGeneratorFigure, 'pointer', oldpointer)
ObsLow=min(cellfun(@(x) min(x(1,1)),SynDatabaseShort));
ObsHigh=max(cellfun(@(x) max(x(1,1)),SynDatabaseShort));
ObsLowest=min(cellfun(@(x) min(x(:,1)),SynDatabaseShort));
ObsHighest=max(cellfun(@(x) max(x(:,1)),SynDatabaseShort));
ht=msgbox(sprintf('Data Generated Successfully! \\ Init Range: [%f, %f]|| Range: [%f, %f]',ObsLow,ObsHigh,ObsLowest,ObsHighest),'Success!');


% AverageY=SumY./DatasetSize;
% for i=1:length(SumY(1,:))
%     SumyOriginal(:,i)=RevSpaceTransform(2,SumY);
% end