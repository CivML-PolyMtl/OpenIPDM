function [pred_vel, pred_vel_var] = TAGI_predict(AnnModel, x_test)
    % Test net
    netT              = AnnModel.net;
    netT.trainMode    = 0;
    netT.batchSize    = 1;
    netT.repBatchSize = 1;
    [netT, statesT, maxIdxT] = network.initialization(netT); 
    normStatT = tagi.createInitNormStat(netT);

    [~, ~, mytest, Sytest] = network.regression(netT, AnnModel.theta_nn, ...
    normStatT, statesT, maxIdxT, x_test, []);

    if AnnModel.net.learnSv
        [mv2a, ~, ~] = act.expFun(mytest(:,2), Sytest(:,2), AnnModel.net.gpu);
        Sytest(:,1) = Sytest(:,1) +  mv2a;
    end
    
    if AnnModel.standardized
        [pred_vel, pred_vel_var] = TAGI_scaler.inverse_transform_norm(mytest(:,1), sqrt(Sytest(:,1)), AnnModel.y_mean_train, AnnModel.y_std_train);
    elseif AnnModel.min_max
        [pred_vel, pred_vel_var] = TAGI_scaler.inverse_transform_min_max_with_var(mytest(:,1), Sytest(:,1),...
                                    AnnModel.y_min, AnnModel.y_max, AnnModel.y_min_range, AnnModel.y_max_range, AnnModel.y_var_factor);
    else
        pred_vel = mytest(:,1);
        pred_vel_var = Sytest(:,1);
    end
    pred_vel = max(min(0, pred_vel),-3);
    pred_vel_var = min(max(pred_vel_var,0.05^2),0.5^2);
end