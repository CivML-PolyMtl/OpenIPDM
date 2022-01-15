% Test input
y_Data_avg = [NaN,70,NaN,NaN,70.77,NaN,NaN,NaN,NaN,67,NaN,NaN,68,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_Data_min = [NaN,41,NaN,NaN,37.77,NaN,NaN,NaN,NaN,55,NaN,NaN,56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_Data_max = [NaN,100,NaN,NaN,98.77,NaN,NaN,NaN,NaN,95,NaN,NaN,97,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
y_perfect = gpuArray(reshape(y_Data_max,1,1,length(y_Data_max)));
y_avg = gpuArray(reshape(y_Data_avg,1,1,length(y_Data_avg)));
y_poor = gpuArray(reshape(y_Data_min,1,1,length(y_Data_min)));
Re = [NaN,1,NaN,NaN,3.87,NaN,NaN,NaN,NaN,48.71,NaN,NaN,4.81,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];
InpecBiase = zeros(1,1,length(Re));
RU = zeros(1,1,length(Re));
ObsYears = [0;1;0;0;1;0;0;0;0;1;0;0;1;0;0;0;0;0;0;0;0;0];
CurrentInspectorObs = zeros(1,length(Re));
A = gpuArray([1,1,0.5;0,1,1;0,0,1]);
C = gpuArray([1,0,0]);
Q = gpuArray([5.0625e-06,1.0125e-05,1.0125e-05;1.0125e-05,2.025e-05,2.025e-05;1.0125e-05,2.025e-05,2.025e-05]);
R = gpuArray(0);
Re = gpuArray(reshape(Re,1,1,length(Re)));
RU = gpuArray(reshape(RU,1,1,length(RU)));
InpecBiase = gpuArray(reshape(InpecBiase,1,1,length(InpecBiase)));
InspBU = InpecBiase;
ObsYears = gpuArray(reshape(ObsYears,1,1,length(ObsYears)));
Ncurve = 4;
init_x_perfect = gpuArray([98.7;-0.8;0]); 
init_x_avg = gpuArray([68;-0.95;0]);
init_x_poor = gpuArray([40;-0.75;0]);
init_V = gpuArray(reshape([9,0,0;0,0.01,0;0,0,0.00249930960072019],3,3,1));
OptBoundsData = [0.005,0.001,0.01;3,0.999,10;3,0.999,10;0.02,0,0.05;0.01,0.001,0.05;0.05,0,0.15];
GlobalCondData = [1,0,25,100;1,1,-50,0;1,0,-25,0];
GPUCompute = 1;
Nsigma = 1;
%% Test #1: Average-health-state-real
% Test 1 verifies the performance of the deterioration model using data for a
% structural element in an average deterioration state
[x, Var, ~, loglik,~,~] = kalman_filter(y_avg, A, C, Q, R, Re...
    , init_x_avg, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
    ,Ncurve,OptBoundsData,GlobalCondData(1,1),GlobalCondData,GPUCompute,Nsigma);
assert(ceil(loglik) == -12);
assert(ceil(x(1,1)) == 68);
assert(ceil(x(2,1)) == 0);
assert(ceil(x(3,1)) == 0);
assert(ceil(x(1,end)) == 61);
assert(ceil(x(2,end)) == 0);
assert(ceil(x(3,end)) == 1);
assert(ceil(Var(1,1,1,1)) == 9);
assert(ceil(Var(2,2,1,1)) == 1);
assert(ceil(Var(3,3,1,1)) == 1);
assert(ceil(Var(1,1,1,end)) == 3);
assert(ceil(Var(2,2,1,end)) == 1);
assert(ceil(Var(3,3,1,end)) == 1);


%% Test #2: Minimum-health-state-real
% Test 2 verifies the performance of the deterioration model using data for a
% structural element in an minimum deterioration state
[x, Var, ~, loglik,~,~] = kalman_filter(y_poor, A, C, Q, R, Re...
    , init_x_poor, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
    ,Ncurve,OptBoundsData,GlobalCondData(1,1),GlobalCondData,GPUCompute,Nsigma);
assert(ceil(loglik) == -38);
assert(ceil(x(1,1)) == 40);
assert(ceil(x(2,1)) == 0);
assert(ceil(x(3,1)) == 0);
assert(ceil(x(1,end)) == 40);
assert(ceil(x(2,end)) == 0);
assert(ceil(x(3,end)) == 1);
assert(ceil(Var(1,1,1,1)) == 9);
assert(ceil(Var(2,2,1,1)) == 1);
assert(ceil(Var(3,3,1,1)) == 1);
assert(ceil(Var(1,1,1,end)) == 1);
assert(ceil(Var(2,2,1,end)) == 1);
assert(ceil(Var(3,3,1,end)) == 1);


%% Test #3: Perfect-health-state-real
% Test 3 verifies the performance of the deterioration model using data for a
% structural element in a perfect health state
[x, Var, ~, loglik,~,~] = kalman_filter(y_perfect, A, C, Q, R, Re...
    , init_x_perfect, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
    ,Ncurve,OptBoundsData,GlobalCondData(1,1),GlobalCondData,GPUCompute,Nsigma);
assert(ceil(loglik) == -10);
assert(ceil(x(1,1)) == 99);
assert(ceil(x(2,1)) == 0);
assert(ceil(x(3,1)) == 0);
assert(ceil(x(1,end)) == 90);
assert(ceil(x(2,end)) == 0);
assert(ceil(x(3,end)) == 1);
assert(ceil(Var(1,1,1,1)) == 9);
assert(ceil(Var(2,2,1,1)) == 1);
assert(ceil(Var(3,3,1,1)) == 1);
assert(ceil(Var(1,1,1,end)) == 5);
assert(ceil(Var(2,2,1,end)) == 1);
assert(ceil(Var(3,3,1,end)) == 1);