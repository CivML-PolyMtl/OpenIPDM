% Test Input

% Observations for a structural element in an average condition
y_Data_avg = [NaN,70,NaN,NaN,70.77,NaN,NaN,NaN,NaN,67,NaN,NaN,68,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% Observations for a structural element in a poor condition
y_Data_min = [NaN,41,NaN,NaN,37.77,NaN,NaN,NaN,NaN,55,NaN,NaN,56,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% Observations for a structural element in a perfect condition
y_Data_max = [NaN,100,NaN,NaN,98.77,NaN,NaN,NaN,NaN,95,NaN,NaN,97,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% Transition matrix 
A = [1,1,0.5,0,0,0;0,1,1,0,0,0;0,0,1,0,0,0;0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1];

% Observation matrix
C = [1,0,0,0,0,0];

% Covariance matrix for the process error of the deterioration model
Q = [2.25e-06,5.63e-06,7.51e-06,0,0,0;
    5.63e-06,1.5e-05,2.25e-05,0,0,0;
    7.51e-06,2.25e-05,4.5e-05,0,0,0;
    0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0];

% Covariance matrix for the process error of the intervention model
Q_r = {[99.9976741869394,0,0,0,0,0;
    0,0.638823842656509,0,0,0,0;
    0,0,0.00249990708482753,0,0,0;
    0,0,0,99.9976741869394,0,0;
    0,0,0,0,0.638823842656509,0;
    0,0,0,0,0,0.00249990708482753]};

% Inspectors' observation variance
Re = [NaN,1,NaN,NaN,3.87,NaN,NaN,NaN,NaN,48.71,NaN,NaN,4.81,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% Inspectors' observation bias
Be = [NaN,0,NaN,NaN,0,NaN,NaN,NaN,NaN,0,NaN,NaN,0,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN,NaN];

% Prior parameters for the interventions model
PriorParam = [0.00671283937681525,1.00103302879139,3,0.0499900434845854,0.0499930955304849,0.149957225722096];

% Initial expected condition for a structural element in an average condition
init_x_avg = [68;-0.95;0;0;0;0];

% Initial expected condition for a structural element in a poor condition
init_x_min = [40;-0.75;0;0;0;0];

% Initial expected condition for a structural element in a perfect condition
init_x_max = [98.7;-0.8;0;0;0;0];

% Initial covariance
init_V = [9,0,0,0,0,0;
    0,0.01,0,0,0,0;
    0,0,0.00249930960072019,0,0,0;
    0,0,0,1.48700305611708,0.00186195874452117,0.000135002413035785;
    0,0,0,0.00186195874452117,0.00966234896454976,5.45609230886886e-05;
    0,0,0,0.000135002413035785,5.46e-05,7.03e-05];

% Transformation function parameter
Ncurve = 4;

% Constrain the deterioration speed using constrained KF model
ConstrainedKF = 1;

% Intervention vector
InterventionCheck = [];
InterventionVector = [0;0;0;0;0;0;0;0;1;0;0;0;0;0;0;0;0;0;0;0;0;0];

% Effect of interventions
InterventionMu_Network = {[15.9329300976325;0.319346922114593;0.000606043592658457]};
InterventionVar_Network = {[1.48700305611708,0.00186195874452117,0.000135002413035785;
    0.00186195874452117,0.00966234896454976,5.45609230886886e-05;
    0.000135002413035785,5.45609230886886e-05,7.03262351532152e-05]};

%% Test #1: Average-health-state-real
% Test 1 verifies the performance of the deterioration model using data for a
% structural element in an average deterioration state
[Ex, VarKF, loglikelihood, Exsmooth, Vsmooth, ~,...
     ~,~,~,~,~,LifeSpanInt]=KF_KS(y_Data_avg,A,C,Q,Q_r,...
     Re,Be,PriorParam,init_x_avg,init_V,Ncurve,ConstrainedKF,InterventionCheck,...
     InterventionVector,InterventionMu_Network,InterventionVar_Network);
assert(ceil(Ex(1,1)) == 68);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 64);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(floor(Exsmooth(1,1)) == 71);
assert(floor(Exsmooth(2,1)) == -1);
assert(floor(Exsmooth(3,1)) == 0);
assert(floor(Exsmooth(1,end)) == 63);
assert(floor(Exsmooth(2,end)) == -1);
assert(floor(Exsmooth(3,end)) == -1);
assert(ceil(Vsmooth(1,1)) == 1);
assert(ceil(Vsmooth(2,2)) == 1);
assert(ceil(Vsmooth(3,3)) == 1);
%% Test #2: Minimum-health-state-real
% Test 2 verifies the performance of the deterioration model using data for a
% structural element in an minimum deterioration state
[Ex, VarKF, loglikelihood, Exsmooth, Vsmooth, ~,...
     ~,~,~,~,~,LifeSpanInt]=KF_KS(y_Data_min,A,C,Q,Q_r,...
     Re,Be,PriorParam,init_x_min,init_V,Ncurve,ConstrainedKF,InterventionCheck,...
     InterventionVector,InterventionMu_Network,InterventionVar_Network);
assert(ceil(Ex(1,1)) == 40);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 49);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(floor(Exsmooth(1,1)) == 41);
assert(floor(Exsmooth(2,1)) == -1);
assert(floor(Exsmooth(3,1)) == -1);
assert(floor(Exsmooth(1,end)) == 48);
assert(floor(Exsmooth(2,end)) == -1);
assert(floor(Exsmooth(3,end)) == -1);
assert(ceil(Vsmooth(1,1)) == 1);
assert(ceil(Vsmooth(2,2)) == 1);
assert(ceil(Vsmooth(3,3)) == 1);
%% Test #3: Perfect-health-state-real
% Test 3 verifies the performance of the deterioration model using data for a
% structural element in a perfect health state
[Ex, VarKF, loglikelihood, Exsmooth, Vsmooth, InterventionMu_Network,...
     InterventionVar_Network,~,~,~,~,LifeSpanInt]=KF_KS(y_Data_max,A,C,Q,Q_r,...
     Re,Be,PriorParam,init_x_max,init_V,Ncurve,ConstrainedKF,InterventionCheck,...
     InterventionVector,InterventionMu_Network,InterventionVar_Network);
assert(ceil(Ex(1,1)) == 99);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 93);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(floor(Exsmooth(1,1)) == 100);
assert(floor(Exsmooth(2,1)) == -1);
assert(floor(Exsmooth(3,1)) == -1);
assert(floor(Exsmooth(1,end)) == 92);
assert(floor(Exsmooth(2,end)) == -1);
assert(floor(Exsmooth(3,end)) == -1);
assert(ceil(Vsmooth(1,1)) == 1);
assert(ceil(Vsmooth(2,2)) == 1);
assert(ceil(Vsmooth(3,3)) == 1);