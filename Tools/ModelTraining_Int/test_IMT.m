% Test input

% Observations for a structural element in a state: min-avg-max
y_min = reshape([NaN;66.8679082110840;NaN;NaN;NaN;NaN;57.2889597661015;62.8664801854814;NaN;NaN;NaN;NaN;NaN;NaN],1,1,14);
y_avg = reshape([NaN;66.8679082110840;NaN;NaN;NaN;NaN;72.3;74;NaN;NaN;NaN;NaN;NaN;NaN],1,1,14);
y_max = reshape([NaN;66.8679082110840;NaN;NaN;NaN;NaN;82;84.86;NaN;NaN;NaN;NaN;NaN;NaN],1,1,14);

% Transition matrix 
A = [1,1,0.5,0,0,0;
    0,1,1,0,0,0;
    0,0,1,0,0,0;
    0,0,0,1,0,0;
    0,0,0,0,1,0;
    0,0,0,0,0,1];

% Inspectors' observation variance: min-avg-max
R_min = reshape([0;17.9167625953763;0;0;0;0;15.9251579974947;19.6936733405843;0;0;0;0;0;0],1,1,14);
R_avg = reshape([0;17.9167625953763;0;0;0;0;15.9251579974947;19.6936733405843;0;0;0;0;0;0],1,1,14);
R_max = reshape([0;17.9167625953763;0;0;0;0;15.9251579974947;19.6936733405843;0;0;0;0;0;0],1,1,14);

% Inspectors' observation bias
Bias = reshape([0;0;0;0;0;0;0;0;0;0;0;0;0;0],1,1,14);

% Transformation function parameter
Ncurve = 4;

% SSM model parameters
param = [0.00415315612822686,1.00007705091361,3,0.0499973324733382,0.0499976538445729,0.149985339623973];
Q = [8.62435291271417e-07,2.15608822817854e-06,2.87478430423806e-06,0,0,0;
    2.15608822817854e-06,5.74956860847611e-06,8.62435291271417e-06,0,0,0;
    2.87478430423806e-06,8.62435291271417e-06,1.72487058254283e-05,0,0,0;
    0,0,0,0,0,0;0,0,0,0,0,0;0,0,0,0,0,0];

% Observation matrix
F =[1 zeros(1,5)];

% Effect of interventions: min-avg-max
InterventionMu_Net_min = [0.210386268117897;0.204895497029854;0.000228625312018737];
InterventionVar_Net_min = [0.0702324899281135,0.000435221277135091,4.22824727048289e-05;
    0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];

InterventionMu_Net_avg = [7;0.204895497029854;0.000228625312018737];
InterventionVar_Net_avg = [3,0.000435221277135091,4.22824727048289e-05;
    0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];

InterventionMu_Net_max = [14;0.204895497029854;0.000228625312018737];
InterventionVar_Net_max = [5,0.000435221277135091,4.22824727048289e-05;
    0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];


% Intervention vector: min-avg-max
InterventionVector_min = reshape([0;0;0;0;0;1;0;0;0;0;0;0;0;0],1,1,14);
InterventionVector_avg = reshape([0;0;0;0;0;1;0;0;0;0;0;0;0;0],1,1,14);
InterventionVector_max = reshape([0;0;0;0;0;1;0;0;0;0;0;0;0;0],1,1,14);
InterventionCheck = 1;

% Constrain the deterioration speed using constrained KF model
ConstrainedKF = 1;

% Covariance matrix for the process error of the intervention model
Q_r_min = [2.41721156793708,0,0,0,0,0;
    0,0.00550982875104675,0,0,0,0;0,0,9.75632695570577e-05,0,0,0;
    0,0,0,2.41721156793708,0,0;0,0,0,0,0.00550982875104675,0;
    0,0,0,0,0,9.75632695570577e-05];
Q_r_avg = [2.41721156793708,0,0,0,0,0;
    0,0.00550982875104675,0,0,0,0;0,0,9.75632695570577e-05,0,0,0;
    0,0,0,2.41721156793708,0,0;0,0,0,0,0.00550982875104675,0;
    0,0,0,0,0,9.75632695570577e-05];
Q_r_max = [2.41721156793708,0,0,0,0,0;
    0,0.00550982875104675,0,0,0,0;0,0,9.75632695570577e-05,0,0,0;
    0,0,0,2.41721156793708,0,0;0,0,0,0,0.00550982875104675,0;
    0,0,0,0,0,9.75632695570577e-05];


% Initial state: min-avg-max (effects)
init_V_min = [17.9167625953763,0,0,0,0,0;
    0,0.0289712926035272,0,0,0,0;0,0,0.00249976538996173,0,0,0;
    0,0,0,0.0702324899281135,0.000435221277135091,4.22824727048289e-05;
    0,0,0,0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    0,0,0,4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];
init_x_min = [66.8679082110840;-0.226123121686847;0;0;0;0];

init_V_avg = [17.9167625953763,0,0,0,0,0;
    0,0.0289712926035272,0,0,0,0;0,0,0.00249976538996173,0,0,0;
    0,0,0,0.0702324899281135,0.000435221277135091,4.22824727048289e-05;
    0,0,0,0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    0,0,0,4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];
init_x_avg = [66.8679082110840;-0.226123121686847;0;0;0;0];

init_V_max = [17.9167625953763,0,0,0,0,0;
    0,0.0289712926035272,0,0,0,0;0,0,0.00249976538996173,0,0,0;
    0,0,0,0.0702324899281135,0.000435221277135091,4.22824727048289e-05;
    0,0,0,0.000435221277135091,0.000290261207466606,1.86234168407862e-05;
    0,0,0,4.22824727048288e-05,1.86234168407862e-05,1.53431404931548e-05];
init_x_max = [66.8679082110840;-0.226123121686847;0;0;0;0];


%% Test #1: Min-health-improvement
% Test 1 verifies the improvement in the health of structural element
[Ex,VarKF,~,loglikely,~,~] = KalmanFilterFun(y_min, A, F, Q, Q_r_min, R_min, Bias, init_x_min, init_V_min,...
    ConstrainedKF,InterventionCheck,InterventionVector_min,...
    InterventionMu_Net_min,InterventionVar_Net_min,param);
assert(ceil(loglikely) == -9);
assert(ceil(Ex(1,1)) == 67);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 61);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(ceil(VarKF(1,1,end)) == 6);
assert(ceil(VarKF(2,2,end)) == 1);
assert(ceil(VarKF(3,3,end)) == 1);

%% Test #2: Avg-health-improvement
% Test 2 verifies the improvement in the health of structural element
[Ex,VarKF,~,loglikely,~,~] = KalmanFilterFun(y_avg, A, F, Q, Q_r_avg, R_avg, Bias, init_x_avg, init_V_avg,...
    ConstrainedKF,InterventionCheck,InterventionVector_avg,...
    InterventionMu_Net_avg,InterventionVar_Net_avg,param);
assert(ceil(loglikely) == -7);
assert(ceil(Ex(1,1)) == 67);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 71);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(ceil(VarKF(1,1,end)) == 7);
assert(ceil(VarKF(2,2,end)) == 1);
assert(ceil(VarKF(3,3,end)) == 1);

%% Test #3: Max-health-improvement
% Test 3 verifies the improvement in the health of structural element
[Ex,VarKF,~,loglikely,~,~] = KalmanFilterFun(y_max, A, F, Q, Q_r_max, R_max, Bias, init_x_max, init_V_max,...
    ConstrainedKF,InterventionCheck,InterventionVector_max,...
    InterventionMu_Net_max,InterventionVar_Net_max,param);
assert(ceil(loglikely) == -8);
assert(ceil(Ex(1,1)) == 67);
assert(ceil(Ex(2,1)) == 0);
assert(ceil(Ex(3,1)) == 0);
assert(ceil(Ex(1,end)) == 80);
assert(ceil(Ex(2,end)) == 0);
assert(ceil(Ex(3,end)) == 0);
assert(ceil(VarKF(1,1,end)) == 7);
assert(ceil(VarKF(2,2,end)) == 1);
assert(ceil(VarKF(3,3,end)) == 1);