[x, Var, ~, loglik,~,~] = kalman_filter(y, A, C, Q, R, Re...
                , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
                ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
TrainSmootherRun();