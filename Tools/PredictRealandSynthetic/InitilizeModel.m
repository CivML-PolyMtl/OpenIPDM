%% Initilize Model
args=varargin;
if args{4} == 0
    is_synthetic = true;
else
    is_synthetic = false;
end
if ~isempty(args)
    GraphMatrix=args{1};                                                    % check if the graphs matrix is provided (Synthetic Analyses)
    FigureID=0;
    FullRun=0;
else
    GraphMatrix=zeros(5);                                                   % default graph matrix
    FullRun=0;
end

% space transformtion
if Pn==1
    MapAxes=2;
    MaxExtendedAxes=37.5;
elseif Pn==2
    MapAxes=1.5;
    MaxExtendedAxes=18.75;
elseif Pn==3
    MapAxes=1.2;
    MaxExtendedAxes=7.5;
elseif Pn==4
    MapAxes=1.1;
    MaxExtendedAxes=3.75;
elseif Pn==5
    MapAxes=1.05;
    MaxExtendedAxes=1.875;
else
    MapAxes=1;
    MaxExtendedAxes=0;
end

% Estimated and True states analysis
ColorLastObs=[];                                                            % For validation with real data (Two databases)
Nsigma=0;                                                                   % The shift in the updated expected value after applying the constraints
y=[NaN y];                                                                  % Add a NaN to the observation so KF start at t=1
ObsYears=[0; ObsYears];                                                     % Make the length of the years vector the same length as the observation vector
if ~isempty(args) && is_synthetic
    FigureID=args{3};
    FullRun=0;
elseif ~isempty(args) && ~is_synthetic
    IndLastObs=find(~isnan(y));                                             % identify the index of the last observation (Real data)
    ObsYears(IndLastObs(end))=0;                                            % Exclude the last observation from the update step in the model
    GraphMatrix=zeros(5);                                                   % Reset the graph matrix to zeros
    ColorLastObs=1;                                                         % trigger the color for the last observation
    FigureID=0;
    FullRun=0;
% elseif ~isempty(args) && length(args)==4
%     if args{4}
%         FullRun=1;
%         GraphMatrix=zeros(5);
%     end
end

InpecBiase=[0 InpecBiase];                                                  % Make the length of the Biase vector the same length as the observation vector
InpecBiaseTrue = [0 InpecBiaseTrue];
OptmInsp=[0 OptmInsp];                                                      % Make the length of the Inspector vector the same length as the observation vector
Re=[0 Re];                                                                  % Make the length of the uncertainnty vector the same length as the observation vector
RU=zeros(1,length(Re));                                                     % Make the length of the current Inspector vector the same length as the observation vector
RunSmoother=find(isnan(y)==0);                                              % find the indexes of observations in the observation vector
RunSmoother=RunSmoother(end);                                               % identify the index of the last observation
if ~isempty(ReTrue)                                                         % check if true data is provided (synthetic database)
    ReTrue=[0 ReTrue];                                                      % Make the length of the true unncertainty vector the same length as the observation vector
    init_xTrue(:,1)=SynDatabaseState{ElementInd}(:,1+yearly(1)-1964);       % get the initial true expected state
    init_VTrue=[1^2,0,0;0,0.1^2,0;0,0,0.1^2];                               % initial uncertainty for the intial state
else
    [~,y]=SpaceTransformation(Pn,y,100.0001,25);                            % transform the observations into the transformed space (Real database)
end

ObsIndexes=find(~isnan(y));                                                 % Identify the indecies of observations (to take the average)
if length(ObsIndexes)==2
    ObsIndexes=ObsIndexes(1:2);
elseif length(ObsIndexes)>2
    ObsIndexes=ObsIndexes(1:3);
else
    ObsIndexes=ObsIndexes(1);
end

% init_x(1)=mean(y(ObsIndexes));                                              % initial condition for the initial state
init_x(1)=y(2);
ConditionVal=mean(y(ObsIndexes));                                                     % initial condition (transformed space)
MaxCondition=100;                                                           % max possible condition (original space)
[Mtrv]=RevSpaceTransform(Pn,ConditionVal,100,25);                                  % initial condition (original space)
DifferenceObs=MaxCondition-Mtrv;                                            % difference between max. condition and initial condition
if is_synthetic
    init_x(2)=TableOfParameters{2,5}*DifferenceObs;                         % expected speed (t=0)
    init_x(3)=TableOfParameters{2,6}*init_x(2);                             % expected acceleration (t=0)
    init_V(1,1,:)=max(TableOfParameters{2,2}.^2,Re(2));                     % condition variance (t=0)
    init_V(2,2,:)=(TableOfParameters{2,3}).^2.*(DifferenceObs)+...          % speed variance (t=0)
        (TableOfParameters{2,7}).^2;
    init_V(3,3,:)=(TableOfParameters{2,4}.^2);                              % acceleration variance (t=0)
else
    init_x(2)=TableOfParameters(1,5)*DifferenceObs;                         % expected speed (t=0) (Real database)
    init_x(3)=TableOfParameters(1,6)*init_x(2);                             % expected acceleration (t=0) (Real database)
    init_V(1,1,:)=max(TableOfParameters(1,2).^2,Re(2));                     % condition variance (t=0) (Real database)
    init_V(2,2,:)=(TableOfParameters(1,3)).^2.*(DifferenceObs)+...          % speed variance (t=0) (Real database)
        (TableOfParameters(1,7)).^2;
    init_V(3,3,:)=(TableOfParameters(1,4).^2);                              % acceleration variance (t=0) (Real database)
end

%% Regression Model
if ~isempty(RegressionModel)
    if kernel_regression
        Kr=1;    
        InirilizedEx=RegressionModel.InirilizedEx;
        InirilizedVar=RegressionModel.InirilizedVar;
        Sigma_W0=RegressionModel.Sigma_W0;
        KernelType=string(RegressionModel.KernelType);
        X_ControlPoints=RegressionModel.X_ControlPoints;
        Kernel_l=RegressionModel.Kernel_l;
        for i=1:length(KernelType)
            Krv=KernelFun(AllAtt(:,i),X_ControlPoints,Kernel_l(i),KernelType(i));
            Kr=Kr.*Krv;
        end
        % Kernel Function
        AKr=Kr./sum(Kr,2);
        init_x(1)=y(2);
        init_x(2,1)=AKr*InirilizedEx;
        init_V(2,2,:)=(AKr*InirilizedVar*AKr' + Sigma_W0^2);
        if init_x(2,1)>0
            init_x(2,1)=0;
        end
    else
        if is_synthetic
            if AnnModel.x_min_maxed
                str_atts = TAGI_scaler.transform_min_max(AllAtt, AnnModel.x_min, AnnModel.x_max, AnnModel.x_min_range , AnnModel.x_max_range);
            elseif AnnModel.x_standardized
                str_atts = TAGI_scaler.transform_norm(AllAtt, AnnModel.x_mean, AnnModel.x_std);
            end
        else
            [x_cat, x_tr] = TAGI_util.split_cat_and_cont_inp(AllAtt, AnnModel.categ_st_att_ind);
            x_cat = TAGI_util.one_hot_encode_inp(x_cat, AnnModel.input_categories);
            % Pre-process inputs
            if AnnModel.x_min_maxed
                x_tr = TAGI_scaler.transform_min_max(x_tr, AnnModel.x_min, AnnModel.x_max, AnnModel.x_min_range , AnnModel.x_max_range);
            elseif AnnModel.x_standardized
                x_tr = TAGI_scaler.transform_norm(x_tr, AnnModel.x_mean, AnnModel.x_std);
            end
            % concatenate the categorical and scaled continuous struct. atrributes
            str_atts = [x_cat, x_tr];
        end

        [pred_vel, pred_vel_var] = TAGI_predict(AnnModel, str_atts);

        if any(isnan(pred_vel)) || any(isnan(pred_vel_var))
            disp('warn! nan values in pred_vel or pred_vel_var')
        end
        init_x(2,1) = pred_vel;
        init_V(2,2,:) = pred_vel_var;
    end
end

dfct_n=@(xts,nts)(nts*exp(-xts.^nts))/gamma(1/nts);                         % dervative of transformation function
d2fct_n=@(xts,nts)-(nts^2*xts^(nts - 1)*exp(-xts^nts))/gamma(1/nts);
load('GlobalCondData');
load('OptBoundsData');
[~, Tdb] = size(y);