StallInit1 = StallVal1;
StallInit2 = StallVal2;
%% Online Inference


%gather data
Y_real = gather(reshape(MdataEngy.YS,[size(MdataEngy.YS,2),size(MdataEngy.YS,3)]));
Inspectorlabel = gather(reshape(MdataEngy.InspectorLabelS,[size(MdataEngy.InspectorLabelS,2),size(MdataEngy.InspectorLabelS,3)]));
Inspector_index = MdataEngy.AllInspectors;

%% lenghts
E = size(Y_real,1);
n   = size(Y_real,2);
n_x = size(3,1);
I = size(Inspector_index,2);

%% initial hidden states

init_x=zeros(3,size(Y_real,1));
init_V = zeros(3,3,length(MdataEngy.init_x));

% acceleration initialisation
init_x(3,:)=0;
init_V(3,3,:)=PARAM(5).^2;

% condition initialisation with max values off 3 first observations
% with mean error ratio
for i=1:size(Y_real,1)
    IndexInit=find(~isnan(gather(Y_real(i,:))));
    if length(IndexInit)==2
        IndexInit=IndexInit(1:2);
    elseif length(IndexInit)>2
        IndexInit=IndexInit(1:3);
    else
        IndexInit=IndexInit(1);
    end
    [valmax, loc]=max(gather(Y_real(i,IndexInit)));
    [tf, loc] = ismember(Inspectorlabel(i,IndexInit(loc)),UpdatedInspectorsData{1}(:,1));
    init_x(1,i) = (1-abs(mean(UpdatedInspectorsData{1}(:,2)))/75)*valmax;  %0.98*valmax;% + 0.4*MdataEngy.init_x(i);%(valmax + valmin) /2 ;% - UpdatedInspectorsData{1}(loc,2);
    init_V(1,1,i)= max(PARAM(3).^2,(UpdatedInspectorsData{1}(loc,3))^2);
end


if isempty(RegressionModel)
    % if regression model is empty
    % speed initialisation
    [~,MAxCondition]=SpaceTransformation(NTr,100,100,25);
    NormObsValues=init_x(1,:);
    NormObsValues(find(NormObsValues>MAxCondition))=MAxCondition;
    DifferenceObs=MAxCondition(end)-NormObsValues;
    init_x(2,:) = 0;
    init_V(2,2,:) = (PARAM(4).^2).*(DifferenceObs)+(PARAM(6).^2);
    
    if OperationIndex>=3
        %with biases
        init_estim = [ 0, 1, PARAM(3), 12]; %[mu_v0, sd(mu_v0), sigma_v0, sd(sigma_v0)]
        [sd, mu_v] = OnlineInference_One(Y_real, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM,init_estim,NTr);
        % saving new inspector's parameters
        UpdatedInspectorsData{1}(:,2) = mu_v';
        UpdatedInspectorsData{1}(:,3) = sd';
        
    else
        % without biases
        init_estim = [ PARAM(3), 12];%[sigma_v0, sd(sigma_v0)]
        [sd ] = OnlineInference_sigma(Y_real, Inspectorlabel, init_x, init_V, A,Q, UpdatedInspectorsData,PARAM,init_estim,MdataEngy);
        % saving new inspector's parameters
        UpdatedInspectorsData{1,1}(:,end) = sd';
    end
    
else
    % speed initialisation with KR parameters
    [AKr,X_ControlPoints]=KR_Prep(MdataEngy.StrucAtt,RegressionModel.KernelType,RegressionModel.Kernel_l,length(RegressionModel.X_ControlPoints)^(1/length(RegressionModel.Kernel_l)));
    init_x(2,:) = AKr*RegressionModel.InirilizedEx;
    init_V(2,2,:)= diag(AKr*RegressionModel.InirilizedVar*AKr') + (RegressionModel.Sigma_W0.^2*ones(size(AKr,1),1));
    
    if OperationIndex>=3
        %with biases
        init_estim = [ 0, 1, PARAM(3), 12]; %  [mu_v0, mu_sd_V0, sd_v0, sd_sd_v0];
        [sd, mu_v] = OnlineInference_short(Y_real, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM,init_estim);
        UpdatedInspectorsData{1}(:,2) = mu_v';
        UpdatedInspectorsData{1}(:,3) = sd';
        
    else
        %without biases
        init_estim = [ PARAM(3), 12]; %  [sd_v0, sd_sd_v0];
        [sd] = OnlineInference_sigma(Y_real, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM,init_estim,MdataEngy);
        UpdatedInspectorsData{1}(:,3) = sd';
    end
end