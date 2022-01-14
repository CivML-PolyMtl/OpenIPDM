% Test Input:
y_perfect = [100,NaN,100,100,NaN,100,NaN,100,NaN,100,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_avg = [75,NaN,77,77,NaN,68,NaN,70,NaN,70,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_poor = [40,NaN,38,39,NaN,32,NaN,35,NaN,33,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
A = [1,1,0.5;0,1,1;0,0,1];
C = [1,0,0];
Q = [5.0625e-06,1.0125e-05,1.0125e-05;1.0125e-05,2.025e-05,2.025e-05;1.0125e-05,2.025e-05,2.025e-05];
R = 0;
Re = [5.86,0,5.94,5.79,0,5.79,0,5.79,0,3.03,0,0,0,0,0,0,0,0];
init_x = [0;0;0];
init_V = [0,0,0;0,0,0;0,0,0];
InpecBiase = [-0.28,0,-0.38,-0.91,0,-0.91,0,-0.91,0,-0.91,0,0,0,0,0,0,0,0];
OptmInsp = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
RU = 1;
InspBU = 0;
ObsYears = [1;0;1;1;0;1;0;1;0;1;0;0;0;0;0;0;0;0];
yearly = [2011,2012,2013,2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024,2025,2026,2027,2028];
InspectorLabel = [5033851,NaN,5038915,5638405,NaN,5638405,NaN,5638405,NaN,6241477,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
StructureInd = 2;
ElementInd = 1;
TableEstimatedParameters = [0.0045,1.02,0.05,0.0245,0,0,0.1499,4];
RegressionModel = [];
AllAtt = [5,65,48.643155,100];
app.BatchMode = 0;
%% Test #1: Average-health-state-real
% Test 1 verifies the performance of the deterioration model using data for a
% structural element in an average deterioration state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,~]=KFsKF(app,y_avg, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -14);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 3);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 76);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 67);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);

%% Test #2: Minimum-health-state-real
% Test 2 verifies the performance of the deterioration model using data for a
% structural element in an minimum deterioration state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,Vvalues]=KFsKF(app,y_poor, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -11);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 3);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 38);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 28);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);

%% Test #3: Perfect-health-state-real
% Test 3 verifies the performance of the deterioration model using data for a
% structural element in a perfect health state
[loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,Vvalues]=KFsKF(app,y_perfect, A, C, Q, R, Re,...
 init_x, init_V,InpecBiase,[],OptmInsp,RU,InspBU,ObsYears,yearly,...
 InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
 [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
assert(ceil(loglik) == -10);
assert(ceil(s_Xsmooth(1,1)) == 2);
assert(ceil(s_Xsmooth(2,1)) == 1);
assert(ceil(s_Xsmooth(3,1)) == 1);
assert(ceil(s_Xsmooth(1,end)) == 2);
assert(ceil(s_Xsmooth(2,end)) == 1);
assert(ceil(s_Xsmooth(3,end)) == 1);
assert(floor(x(1,1)) == 106);
assert(floor(x(2,1)) == -1);
assert(floor(x(3,1)) == -1);
assert(floor(x(1,end)) == 100);
assert(floor(x(2,end)) == -1);
assert(floor(x(3,end)) == -1);