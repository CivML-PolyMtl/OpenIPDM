clc
dt = 1.5;
A=[1 dt dt^2/2;0 1 dt;0 0 1];                                              % transition for the kinimatic model                                                     % transition matrix (with intervention)
c =[1, 0, 0]; 
constrain_vector = [1, 0, 0];
A = constr_component_a(A,constrain_vector,c);

mu = [0.1; 0.2; 0.3];
covariance = [0.1, 0.2, 0.3; 
              0.2, 0.4, 0.5; 
              0.3, 0.5, 0.8];

covariance = covariance * covariance';

lower_bounds = [-0.1, -inf, -inf];
upper_bounds = [0.1, inf, inf];

[x, p] = relu_constrain_step(mu, covariance, constrain_vector, lower_bounds, upper_bounds, false);
disp(x)
disp(p)

