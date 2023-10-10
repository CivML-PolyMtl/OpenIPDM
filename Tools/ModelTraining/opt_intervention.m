function [tot_ll,InterventionMu_Network,InterventionVar_Network]= ...
    opt_intervention(Cindex,InspectorsID,InspectorsData,Re,Be,y,...
    param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
    InterventionVector,model_i,AnnModel,AllAtt,int_type)
tot_ll=0;
% Prior Knowledge
if int_type == 3
    InterventionMu_Network=[10 0.4 0];
elseif int_type == 2
    InterventionMu_Network=[5 0.4 0];
else
    InterventionMu_Network=[0 0.8 0];
end
InterventionVar_Network=diag([ModelParamLocal(4)^2 ModelParamLocal(5)^2 ModelParamLocal(6)^2]);%[5^2 0.3^2 0.05^2]);%

Q_r=diag([ModelParamLocal(1)^2 ModelParamLocal(2)^2 ModelParamLocal(3)^2]);%diag([2^2 0.1^2 0.01^2]);
Q_r=blkdiag(Q_r,Q_r);

for i=1:length(Cindex)
    for j=1:length(InspectorsID(1,i,:))
        InspID=find(InspectorsID(1,Cindex(i),j)==InspectorsData{1}(:,1));
        if ~isempty(InspID)
            Re(1,Cindex(i),j)=InspectorsData{1}(InspID,3)^2;
            Be(1,Cindex(i),j)=InspectorsData{1}(InspID,2);
        end
    end
    % Initial values
    init_x(1)=y(1,Cindex(i),2);
    init_x(3:6,1)=zeros(4,1);
    init_V(1,1)=max(param(3).^2,Re(1,Cindex(i),2));
    init_V(3,3)=param(5).^2;
    init_V(4:6,4:6)=InterventionVar_Network;
    % BNN Regression
    
    if model_i==2
            if AnnModel.categ_st_att_ind > 0
                [x_cat, x_tr] = TAGI_util.split_cat_and_cont_inp(AllAtt(i,:), AnnModel.categ_st_att_ind);
                % one-hot encode
                x_cat = TAGI_util.one_hot_encode_inp(x_cat,AnnModel.input_categories);
            else
                x_tr = AllAtt(i,:);
            end
        
        % x_tr = StrucAtt_cpu;
        
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
        init_x(2,:) = pred_vel;
        init_V(2,2) = pred_vel_var;
        if isnan(init_x(2)) || isnan(init_V(2,2))
            message = 'Model training has failed due to KR model providing a NaN deterioration speed estimate. Consider adjusting the KR framework input and/or parameters.';
            msgbox(message, 'Model training is incomplete')
        end
    else
        init_x(2,:)=0;
        [MAxCondition]=RevSpaceTransform(Ncurve,100,100.001,25);
        DifferenceObs=MAxCondition-init_x(1);
        init_V(2,2)=(param(4).^2).*(DifferenceObs)+(param(6).^2);
    end
    try 
    [~, ~, loglikelihood, ~, ~, InterventionMu_Network,...
        InterventionVar_Network]=HandCoded_KF_Network(y(1,Cindex(i),:),A,F,Q,Q_r,...
        Re(1,Cindex(i),:),Be(1,Cindex(i),:),param,init_x,init_V,Ncurve,ConstrainedKF,...
        InterventionCheck,InterventionVector(1,Cindex(i),:),...
        InterventionMu_Network,InterventionVar_Network);
    catch
        init_x
    end
    tot_ll= tot_ll + loglikelihood;
end
end