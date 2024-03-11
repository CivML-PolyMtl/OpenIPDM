dt = 1.5;
A = [1, dt, (dt^2)/2; 0, 1, dt; 0, 0, 1];
c = [1, 0, 0];
constrain_vector = [0, 1, 0];

A = constr_component_a(A,constrain_vector,c);

disp(A)

