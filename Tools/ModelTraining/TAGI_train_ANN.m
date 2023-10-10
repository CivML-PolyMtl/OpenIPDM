% Name the experiment
modelName = 'SSM_BNN';
dataName = 'network_level';
AnnModel.dataName = dataName;
AnnModel.modelName = modelName;

global total_time_TAGI;

% Input for TAGI
StrucAtt_train = reshape(gather(StrucAtt),size(StrucAtt,2),size(StrucAtt,3));
StrucAtt_valid = reshape(gather(MdataEngy.ModelValid.StrucAtt),size(MdataEngy.ModelValid.StrucAtt,2),size(MdataEngy.ModelValid.StrucAtt,3));
% [x_train, x_val, AnnModel] = preprocess_tagi_inputs(StrucAtt_train, StrucAtt_valid, AnnModel);
x_train_val = StrucAtt_train;
x_test = StrucAtt_valid;

% one-hot encode
% categories = unique(x_train_val_cat);
categories = MdataEngy.StrucAttCategories;
AnnModel.input_categories = gather(categories);

if AnnModel.categ_st_att_ind > 0
    [x_train_val_cat, x_train_val_cont] = TAGI_util.split_cat_and_cont_inp(x_train_val, AnnModel.categ_st_att_ind);
    % one-hot encode
    x_train_val_cat = TAGI_util.one_hot_encode_inp(x_train_val_cat,categories);
else
    x_train_val_cont = x_train_val;
end

% Pre-process inputs
% x_min_range = -1;
% x_max_range = 1;
AnnModel.x_min_maxed = false;
AnnModel.x_standardized = true;

if AnnModel.x_min_maxed
    [x_train_val_cont, x_min, x_max] = TAGI_scaler.fit_transform_min_max(x_train_val_cont, x_min_range, x_max_range);
    AnnModel.x_min = x_min;
    AnnModel.x_max = x_max;
    AnnModel.x_min_range = x_min_range;
    AnnModel.x_max_range = x_max_range;
elseif AnnModel.x_standardized
    [x_train_val_cont, x_mean, x_std] = TAGI_scaler.fit_transform_norm(x_train_val_cont);
    AnnModel.x_mean = x_mean;
    AnnModel.x_std = x_std;
end

% concatenate the categorical and scaled continuous struct. atrributes
if AnnModel.categ_st_att_ind > 0
    x_train_val_cont = [x_train_val_cat, x_train_val_cont];
end


%% Global validation set preparation 
if AnnModel.categ_st_att_ind > 0 
    [x_test_cat, x_test_cont] = TAGI_util.split_cat_and_cont_inp(x_test, AnnModel.categ_st_att_ind);
    x_test_cat = TAGI_util.one_hot_encode_inp(x_test_cat,AnnModel.input_categories);
else
    x_test_cont = x_test;
end

if AnnModel.x_min_maxed
    x_test_cont = TAGI_scaler.transform_min_max(x_test_cont, x_min, x_max, x_min_range, x_max_range);
elseif AnnModel.x_standardized
    x_test_cont = TAGI_scaler.transform_norm(x_test_cont, x_mean, x_std);
end
if AnnModel.categ_st_att_ind > 0 
    x_test_cont = [x_test_cat, x_test_cont]; % concatenate the categorical and scaled continuous struct. atrributes
end


% % 
% y_min_range = -0.5;
% y_max_range = 0.5;
% AnnModel.y_min_range = y_min_range;
% AnnModel.y_max_range = y_max_range;

%% TAGI set-up
define_tagi_architecture();

%% Regression loop and validation set preparation 
prep_regress_loop_and_val_set();

%% Initialize the state
[init_x, init_V] = define_init_state.plain(y, Re, param);

%% Global loop
time_start = tic;
while RegressLoop<MultiPass && ~Converge
    
    DifferenceObs=MaxCondition(end)-NormObsValues;
    
    fprintf('\nRegress loop: %d\n', RegressLoop);

    % Kalman Filter and Smoother to generate training data for our 'hyper' model
    try
        if RegressLoop == 1 && isempty(AnnModel.theta_nn)
            init_x(2,:)= 0 ;
            init_V(2,2,:)=(3*param(4).^2).*(DifferenceObs)+(param(6).^2);
            
            KF_KS_train();

            InitialCond=gather(ExSmooth(1,:,2));
            % update condition values with smoothed values
            InitialCond(find(InitialCond>MaxCondition(end)))=MaxCondition(end);
            DifferenceObs=MaxCondition(end)-InitialCond;

            init_V(2,2,:)=(3*param(4).^2).*(DifferenceObs)+(param(6).^2);
            
            KF_KS_train();
        else
            [pred_vel, pred_vel_var] = TAGI_predict(AnnModel, x_train_val_cont);
    
            [init_x, init_V] = define_init_state.with_BNN(pred_vel, pred_vel_var, y, Re, param);

            init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);

            KF_KS_train();

            InitialCond=gather(ExSmooth(1,:,2));
            % update condition values with smoothed values
            InitialCond(find(InitialCond>MaxCondition(end)))=MaxCondition(end);
            DifferenceObs=MaxCondition(end)-InitialCond;

            init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);

            KF_KS_train();
        end
    catch
        disp(' ')
        disp('Constraints failiure during training data generation')
        SumLoglikValidation=-99^10;
        loglik=-99^10;
        break;
    end
    % retrieve initial speed mean and variance from KS
    Sy_train_val=reshape(gather(VarSmooth(2,2,:,1)),length(VarSmooth(2,2,:,1)),1);
    y_train_val=gather(ExSmooth(2,:,1)'); % ind = 1: condition || ind = 2: speed || ind = 3 : acc    
    AnnModel.y_tr = y_train_val;

    [y_train_val, y_mean_train, y_std_train] = TAGI_scaler.fit_transform_norm(y_train_val);
    Sy_train_val = Sy_train_val./(y_std_train.^2); % scale the variance used for training as well
    AnnModel.standardized = true;
    AnnModel.y_mean_train = y_mean_train;
    AnnModel.y_std_train = y_std_train;
    AnnModel.min_max = false; 
    
    [x_train, y_train, Sy_train, x_val, y_val] = TAGI_train_test_split(x_train_val_cont, y_train_val, Sy_train_val, split_ratio, seed, shuffle);
    
    % Train
    [theta, LL_val, exceeded_reinit_attempts] = TAGI_train(net, x_train, y_train,  Sy_train, x_val, y_val, seed, AnnModel.theta_nn);
    
    if ~exceeded_reinit_attempts
        AnnModel.net = net;
        AnnModel.theta_nn = theta;

        %% Validation with validation set
        [pred_vel, pred_vel_var] = TAGI_predict(AnnModel, x_test_cont);

        [init_x_valid, init_V_valid] = define_init_state.with_BNN(pred_vel, pred_vel_var, y_valid, Re_Valid, param);

        % Kalman Filter on the validation set
        KF_validation_set()
        if break_
            break
        end

        if ~isempty(AnnModel_loaded) && RegressLoop == 1
            [pred_vel, pred_vel_var] = TAGI_predict(AnnModel_loaded, x_test_cont);
            [init_x_valid, init_V_valid] = define_init_state.with_BNN(pred_vel, pred_vel_var, y_valid, Re_Valid, param);

            % Kalman Filter on the validation set
            KF_validation_set()     
            if break_
                break
            end
        else
            LogLikValidation_loaded = -inf;
        end

        %% Stopping Criteria    
        SumLL(RegressLoop,1)=gather(sum(LogLikValidation(:)));
        fprintf('Regress loop #%d LL validation: %d\n', RegressLoop, SumLL(RegressLoop,1));

        regress_loop_convergence_check();

        RegressLoop=RegressLoop+1;
    else
        Converge = 1;
        Best_RegressLoop = 0;
    end
    %% If converged, evaluate the estimated parameters
    if RegressLoop>=MultiPass || Converge
        report_convergence();
             
        [pred_vel, pred_vel_var]= TAGI_predict(AnnModel, x_train_val_cont);
        
        [init_x, init_V] = define_init_state.with_BNN(pred_vel, pred_vel_var, y, Re, param);
        try
            KF_KS_train()
        catch
           disp('LL eval of train set failed (at the end of regress loop convergence)') 
        end
    end
end

