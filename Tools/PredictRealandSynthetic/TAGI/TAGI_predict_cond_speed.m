function [pred_state, pred_state_var, pred_cov] = TAGI_predict_cond_speed(AnnModel, x_test)
    % Test net
    netT              = AnnModel.net;
    netT.trainMode    = 0;
    netT.batchSize    = 1;
    netT.repBatchSize = 1;
    [netT, statesT, maxIdxT] = network.initialization(netT); 
    normStatT = tagi.createInitNormStat(netT);

    [~, ~, pred_mean, pred_var] = network.regression(netT, AnnModel.theta_nn, ...
    normStatT, statesT, maxIdxT, x_test, []);

    % Transform from the 'Cholesky space' to the 'variance space'
    mv2_test = zeros(size(pred_mean(:,netT.nl+1:end)));
    Sv2_test = zeros(size(pred_mean(:,netT.nl+1:end)));

    for i=1:size(pred_mean,1)
        ma = pred_mean(i,:);
        Sa = pred_var(i,:);
        [~, mLa] = tagi.detachMeanVar(ma, netT.nl,...
            netT.nLchol, netT.batchSize, netT.repBatchSize);
        [~, SLa] = tagi.detachMeanVar(Sa, netT.nl,...
            netT.nLchol, netT.batchSize, netT.repBatchSize);
        % transform the diagonal elements into positive domain
        [mLa_, SLa_, ~] = agvi.transform_chol_vec(mLa, SLa, netT.gpu);
        % retrieve the variance parameters from the cholesky vectors
        [mv2, Sv2] = agvi.chol_to_v2(mLa_, SLa_);
        mv2_test(i,:) = mv2;
        Sv2_test(i,:) = diag(Sv2);
    end

    var_indices = [1, 3];

    for i=1:netT.nl
        pred_var(:,i) = pred_var(:,i) + mv2_test(:,var_indices(i));
    end

    pred_corr =  mv2_test(:,2)./sqrt((mv2_test(:,1)).*sqrt(mv2_test(:,3)));
    
    if AnnModel.standardized
        [pred_state, pred_state_var] = TAGI_scaler.inverse_transform_norm(pred_mean(:,1:2), sqrt(pred_var(:,1:2)), AnnModel.y_mean_train, AnnModel.y_std_train);
        [~,pred_aleatory_var]=TAGI_scaler.inverse_transform_norm(pred_mean(:,1:2), sqrt(mv2_test(:,var_indices)), AnnModel.y_mean_train, AnnModel.y_std_train);
        pred_cov = pred_corr .* sqrt(pred_aleatory_var(:,1)) .* sqrt(pred_aleatory_var(:,2));
    elseif AnnModel.min_max
        [pred_state, pred_state_var] = TAGI_scaler.inverse_transform_min_max_with_var(pred_mean(:,1:2), pred_var(:,1:2),...
                                    AnnModel.y_min, AnnModel.y_max, AnnModel.y_min_range, AnnModel.y_max_range, AnnModel.y_var_factor);
    else
        pred_state = pred_mean(:,1:2);
        pred_state_var = pred_var(:,1:2);
    end
end