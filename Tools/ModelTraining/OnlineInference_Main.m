StallInit1 = StallVal1;
StallInit2 = StallVal2;
%% Online Inference


%gather data
Y_real = gather(reshape(ElementData.YS,[size(ElementData.YS,2),size(ElementData.YS,3)]));
Inspectorlabel = gather(reshape(ElementData.InspectorLabelS,[size(ElementData.InspectorLabelS,2),size(ElementData.InspectorLabelS,3)]));
Inspector_index = ElementData.AllInspectors;

%% lenghts
E = size(Y_real,1);
n   = size(Y_real,2);
n_x = size(3,1);
I = size(Inspector_index,2);

%% initial hidden states
init_x=zeros(3,size(Y_real,1));
init_V = zeros(3,3,length(ElementData.init_x));

% acceleration initialisation
init_x(3,:)=0;
init_V(3,3,:)=PARAM(5).^2;

% condition initialisation with max values off 3 first observations
% with mean error ratio
NormObsValues = [];
for i=1:size(Y_real,1)
    IndexInit=find(~isnan(gather(Y_real(i,:))));
    if length(IndexInit)==2
        IndexInit=IndexInit(1:2);
    elseif length(IndexInit)>2
        IndexInit=IndexInit(1:3);
    else
        IndexInit=IndexInit(1);
    end
    NormObsValues = [NormObsValues, Y_real(i, IndexInit(1))];
    [valmax, loc]=max(gather(Y_real(i,IndexInit)));

    [tf, loc] = ismember(Inspectorlabel(i,IndexInit(loc)),InspectorsData{1}(:,1));
    init_x(1,i) = (1-abs(mean(InspectorsData{1}(:,2)))/75)*valmax;  %0.98*valmax;% + 0.4*MdataEngy.init_x(i);%(valmax + valmin) /2 ;% - UpdatedInspectorsData{1}(loc,2);
    init_V(1,1,i)= max(PARAM(3).^2,(InspectorsData{1}(loc,3))^2);
end

if (isempty(AnnModel) || isempty(AnnModel.theta_nn)) % the or condition is for when the AnnModel has not yet been trained
    % speed initialisation using plain SSM -- i.e. zero speed prior
    [~,MAxCondition]=SpaceTransformationVec(Ncurve,100,100,25);
    %NormObsValues=init_x(1,:); % potential bug
    NormObsValues(find(NormObsValues>MAxCondition))=MAxCondition;
    DifferenceObs = MAxCondition(end)-NormObsValues;
    init_x(2,:) = 0;
    init_V(2,2,:) = (PARAM(4).^2).*(DifferenceObs)+(PARAM(6).^2);
    
    %with biases
    
    [sd, mu_v, cov_params] = OnlineInference_One(Y_real, Inspectorlabel, init_x, init_V, A, Q, InspectorsData,PARAM,init_estim,Ncurve);
    % saving new inspector's parameters
    InspectorsData{1}(:,2) = mu_v';
    InspectorsData{1}(:,3) = sd';
        
else
    str_atts = ElementData.StrucAtt;
    str_atts = reshape(gather(str_atts), size(str_atts,2), size(str_atts,3));
    
    if AnnModel.categ_st_att_ind > 0 
        [x_cat, x_tr] = TAGI_util.split_cat_and_cont_inp(str_atts, AnnModel.categ_st_att_ind);
        x_cat = TAGI_util.one_hot_encode_inp(x_cat,AnnModel.input_categories);
    else
        x_tr = str_atts;
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
    
    [pred_vel, pred_vel_var] = TAGI_predict(AnnModel, x_tr);
    
    init_x(2,:) = pred_vel;
    init_V(2,2,:) = pred_vel_var;

    %with biases
    [sd, mu_v, cov_params] = OnlineInference_short(Y_real, Inspectorlabel, init_x, init_V, A, Q, InspectorsData,PARAM,init_estim);
    InspectorsData{1}(:,2) = mu_v';
    InspectorsData{1}(:,3) = sd';

end