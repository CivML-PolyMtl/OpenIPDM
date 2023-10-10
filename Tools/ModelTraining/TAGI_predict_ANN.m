% Function returns the loglik of KF, using TAGI to initialize speeds

StrucAtt_cpu = reshape(gather(StrucAtt),size(StrucAtt,2),size(StrucAtt,3));

if AnnModel.categ_st_att_ind > 0 
    [x_cat, x_tr] = TAGI_util.split_cat_and_cont_inp(StrucAtt_cpu, AnnModel.categ_st_att_ind);
    x_cat = TAGI_util.one_hot_encode_inp(x_cat,AnnModel.input_categories);
else
    x_tr = StrucAtt_cpu;
end

% Pre-process inputs
if AnnModel.x_min_maxed
    x_tr = TAGI_scaler.transform_min_max(x_tr, AnnModel.x_min, AnnModel.x_max, AnnModel.x_min_range , AnnModel.x_max_range);
elseif AnnModel.x_standardized
    x_tr = TAGI_scaler.transform_norm(x_tr, AnnModel.x_mean, AnnModel.x_std);
end
% concatenate the categorical and scaled continuous struct. atrributes
if AnnModel.categ_st_att_ind > 0 
    x_tr = [x_cat, x_tr];
end

% Use TAGI to get initial speed prediction
[pred_vel, pred_vel_var] = TAGI_predict(AnnModel, x_tr);

[init_x, init_V] = define_init_state.with_BNN(pred_vel, pred_vel_var, y, Re, param);

[~, ~, ~, loglik,~,~] = kalman_filter(y, A, C, Q, R, Re...
    , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
    ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
