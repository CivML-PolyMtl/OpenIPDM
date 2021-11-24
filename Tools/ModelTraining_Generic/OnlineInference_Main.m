StallInit1 = StallVal1;
StallInit2 = StallVal2;
%% Online Inference


%gather data
Y_real = gather(reshape(MdataEngy.YS,[size(MdataEngy.YS,2),size(MdataEngy.YS,3)]));
Inspectorlabel = gather(reshape(MdataEngy.InspectorLabelS,[size(MdataEngy.InspectorLabelS,2),size(MdataEngy.InspectorLabelS,3)]));
Inspector_index = MdataEngy.AllInspectors;

%lenghts
E = size(Y_real,1);
n   = size(Y_real,2);
n_x = 3;
I = size(Inspector_index,2);

%initial hidden states
% 			init_x(1,:)=MdataEngy.init_x;

init_x=zeros(3,size(Y_real,1));
init_V = zeros(3,3,length(MdataEngy.init_x));

init_x(3,:)=0;

init_x(1,:) = Y_real(:,2)';
for e = 1:E
    num_ins   = Inspectorlabel(e,2);
    index_ins = find(UpdatedInspectorsData{1,1}==num_ins);
    sigma_ins = UpdatedInspectorsData{1,1}(index_ins,end);
    init_V(1,1,e)=max(PARAM(2).^2,sigma_ins.^2);
end

init_V(3,3,:)=PARAM(5).^2;

% min ou max
%             for i=1:size(Y_real,1)
%                 IndexInit=find(~isnan(gather(Y_real(i,:))));
%                 if length(IndexInit)==2
%                     IndexInit=IndexInit(1:2);
%                 elseif length(IndexInit)>2
%                     IndexInit=IndexInit(1:3);
%                 else
%                     IndexInit=IndexInit(1);
%                 end
%                 [valmax, loc]=max(gather(Y_real(i,IndexInit)));
%                 [tf, loc] = ismember(Inspectorlabel(i,IndexInit(loc)),UpdatedInspectorsData{1}(:,1));
%                 init_x(1,i) = valmax ;% - UpdatedInspectorsData{1}(loc,2);
%                 init_V(1,1,i)=max(PARAM(3).^2,(UpdatedInspectorsData{1}(loc,3))^2);
%             end


if isempty(RegressionModel)
    % initialisation condition with first value
    %                 init_x(1,:) = MdataEngy.init_x ; % Y_real(:,2)';
    %                 init_x(1,:) = Y_real(:,2)';
    %                 for i=1:size(Y_real,1)
    %                     IndexInit=find(~isnan(gather(Y_real(i,:))));
    %                     if length(IndexInit)==2
    %                         IndexInit=IndexInit(1:2);
    %                     elseif length(IndexInit)>2
    %                         IndexInit=IndexInit(1:3);
    %                     else
    %                         IndexInit=IndexInit(1);
    %                     end
    %                     init_x(1,i)=max(gather(Y_real(i,IndexInit)));
    %                 end
    %                 %initialisation speed condition with p1 and p2
    [~,MAxCondition]=SpaceTransformation(NTr,100,100,25);
    NormObsValues=init_x(1,:);
    NormObsValues(find(NormObsValues>MAxCondition))=MAxCondition;
    DifferenceObs=MAxCondition(end)-NormObsValues; % first time
    init_x(2,:) = 0;
    init_V(2,2,:) = (PARAM(4).^2).*(DifferenceObs)+(PARAM(6).^2);
    
    if OperationIndex>=3
        %with biases
        init_estim = [ 0, 2, PARAM(3), 8];
        % first run
        [sd, mu_v] = OnlineInference_One(Y_real, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM,init_estim,NTr);
        % save
        UpdatedInspectorsData{1}(:,2) = mu_v';
        UpdatedInspectorsData{1}(:,3) = sd';
        
    else
        % without biases
        init_estim = [ PARAM(3), 8];
        % first run
        [sd ] = OnlineInference_sigma(Y_real, Inspectorlabel, init_x, init_V, A,Q, UpdatedInspectorsData,PARAM,init_estim,MdataEngy);
        % save
        UpdatedInspectorsData{1,1}(:,end) = sd';
    end
    
else
    %                 init_x(1,:) = Y_real(:,2)';
    %                 % initialisation condition with max value.
    %                 for i=1:size(Y_real,1)
    %                     IndexInit=find(~isnan(gather(Y_real(i,:))));
    %                     if length(IndexInit)==2
    %                         IndexInit=IndexInit(1:2);
    %                     elseif length(IndexInit)>2
    %                         IndexInit=IndexInit(1:3);
    %                     else
    %                         IndexInit=IndexInit(1);
    %                     end
    %                     init_x(1,i)=max(gather(Y_real(i,IndexInit)));
    % %                 end
    % UpdatedInspectorsData{1,1}(:,3) = PARAM(3);
    % UpdatedInspectorsData{1,1}(:,2) = 0;
    % initialisation with KR parameter estiation for the speed condition.
    [AKr,X_ControlPoints]=KR_Prep(MdataEngy.StrucAtt,RegressionModel.KernelType,RegressionModel.Kernel_l,length(RegressionModel.X_ControlPoints)^(1/length(RegressionModel.Kernel_l)));
    init_x(2,:) = AKr*RegressionModel.InirilizedEx;
    init_V(2,2,:)= diag(AKr*RegressionModel.InirilizedVar*AKr') + (RegressionModel.Sigma_W0.^2*ones(size(AKr,1),1));
    % compute the max of the first 2:3 observations
    
    if OperationIndex>=3
        %with biases
        init_estim = [ 0, 2, PARAM(3), 12]; %  [mu_v0, mu_sd_V0, sd_v0, sd_sd_v0];
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