%% Initialize the state 
init_x(1,:)=gather(y(1,:,2));
if OptProcedure==1
    init_x(2,:)=0;%param(6)*(DifferenceObs);%+param(8);
    init_x(3,:)=0;%param(7)*init_x(2,:);
else
    init_x(2,:)=-0.01;
    init_x(3,:)=-0.001;
end
if OptProcedure==1
    if GPUCompute
        init_V=zeros(3,3,length(IniObsAvg(1,:)),'gpuArray');
        
    else
        init_V=zeros(3,3,1);
    end
    if OptLevel~=1
        param(3)=param(2);
    end
    init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
    init_V(3,3,:)=param(5).^2;%(param(5).^2).*init_V(2,2,:)+init_V(2,2,:).*param(7)^2+ 
end
tic
% Prepare Akr / Control Points
[AKr,X_ControlPoints]=KR_Prep(StrucAtt,KernelType,Kernel_l,Ncp);

% Initial value for Hidden Control Points
if args{6}<=1
    InirilizedEx=ones(Ncp^length(KernelType),1)*0;
end
% Regression Loop (Multipass)
RegressLoop=1;
SumLL=zeros(MultiPass,1);
Converge=0;
TolRegress=10^-3;
TotalTimeSteps=size(y,3);
OriginalValues=[25:1:100];
[~,MAxCondition]=SpaceTransformation(Ncurve,OriginalValues,100,25);
NormObsValues=init_x(1,:);
NormObsValues(find(NormObsValues>MAxCondition(end)))=MAxCondition(end);
% init_x(1,:)=IniObsAvg(1,:);%gather(y(1,:,2));

% Validation Config.
if args{6}<=1
    StrucAtt_Valid=MdataEngy.ModelValid.StrucAtt;
    y_valid=MdataEngy.ModelValid.YS;
    Re_Valid=MdataEngy.ModelValid.ReS;
    R_valid=MdataEngy.ModelValid.RS;
    RU_valid=MdataEngy.ModelValid.RUS;
    InpecBiase_valid=MdataEngy.ModelValid.InpecBiaseS;
    InspBU_valid=MdataEngy.ModelValid.InspBUS;
    ObsYears_valid=MdataEngy.ModelValid.ObsYearsS;
    InspectorsObs_valid=MdataEngy.ModelValid.AllInspectors;
    for i=1:length(InspectorsID)
        ReReshape_valid=reshape(Re_Valid,[size(Re_Valid,2),size(Re_Valid,3)]);
        ReReshape_valid(find(InspectorsObs_valid{i}))=EngBiasData(i,end).^2;
        Re_Valid(1,:,:)=ReReshape_valid;
    end
    CurrentInspectorObs_valid=zeros(size(InspectorsObs_valid{1}),'gpuArray');
    init_V_valid=zeros(3,3,length(MdataEngy.ModelValid.RS),'gpuArray');
    if ~isempty(CurrentInspectorID)
        RUReshape_valid=reshape(Re_Valid,[size(Re_Valid,2),size(Re_Valid,3)]);
        CurrentInspectorIndex_valid=find(CurrentInspectorID==InspectorsID);
        CurrentInspectorObs=InspectorsObs_valid{CurrentInspectorIndex_valid};
        RUReshape_valid(find(CurrentInspectorObs))=(CurrentInspectorParam(1)).^2;
        Re_Valid(1,:,:)=RUReshape_valid;
    end
end
%% Loop

while RegressLoop<MultiPass && ~Converge
%     T_0=cputime;
    
    DifferenceObs=MAxCondition(end)-NormObsValues;
    if args{6}<=1
        % re-initialize the variance
        InirilizedVar=diag(ones(length(InirilizedEx),1).*2.^2);
        init_x(2,:)=AKr*InirilizedEx;
        init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);
        % Kalman Filter and Smoother
        try
        [x, Var, ~, ~,~,~] = kalman_filter(y, A, C, Q, R, Re...
            , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
            ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
        TrainSmootherRun();      
        if GPUCompute
            InitialCond=gather(ExSmooth(1,:,2));
            % update condition values with smoothed values
            InitialCond(find(InitialCond>MAxCondition(end)))=MAxCondition(end);
            DifferenceObs=MAxCondition(end)-InitialCond;
        else
            InitialCond=ExSmooth(1,2);
            InitialCond(find(InitialCond>MAxCondition(end)))=MAxCondition(end);
            DifferenceObs=MAxCondition(end)-InitialCond;
        end
        init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);
        [x, Var, ~, ~,~,~] = kalman_filter(y, A, C, Q, R, Re...
            , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
            ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
        TrainSmootherRun();      
        catch
            disp(' ')
            disp('Constraints Failiur')
            SumLoglikValidation=-99^10;
            loglik=-99^10;
            break;
        end
        % Variance from Kalman Smoother
        VarSmootherVal=diag(reshape(gather(VarSmooth(2,2,:,1)),length(VarSmooth(2,2,:,1)),1));
        XpredSpeed=init_x(2,:)';
        VpredKr=(AKr*InirilizedVar*AKr') + diag(Sigma_W0^2*ones(size(AKr,1),1));
        StoreInirilizedEx=InirilizedEx;
        StoreInirilizedVar=InirilizedVar;
        Xsmooth=ExSmooth;
        % Update hidden control points
        J=InirilizedVar*AKr'/(VpredKr);
        InirilizedEx=InirilizedEx+J*gather(Xsmooth(2,:,1)'-XpredSpeed);
        InirilizedVar=InirilizedVar+J*gather(VarSmootherVal-VpredKr)*J';
        T_T=toc;
        %% Validation with validation set
        % Prepare Akr / Control Points
        [AKr_Valid,~]=KR_Prep(StrucAtt_Valid,KernelType,Kernel_l,Ncp);
        
        % initial speed state
        init_x_valid(1,:)=gather(y_valid(1,:,2));%MdataEngy.ModelValid.init_x(1,:);
        init_x_valid(2,:)=AKr_Valid*InirilizedEx;
        init_x_valid(3,:)=0;
        init_V_valid(1,1,:)=max(param(3).^2,Re_Valid(1,:,2));
        init_V_valid(2,2,:)=diag(AKr_Valid*InirilizedVar*AKr_Valid')' + (Sigma_W0^2*ones(size(AKr_Valid,1),1))';%(param(4).^2).*(DifferenceObs)+(param(6).^2);
        init_V_valid(3,3,:)=param(5).^2;
        
        % Kalman Filter and Smoother
        try
            [x, Var, ~, LogLikValidation,~,~] = kalman_filter(y_valid, A, C, Q, R_valid, Re_Valid...
                , init_x_valid, init_V_valid,InpecBiase_valid,CurrentInspectorObs_valid,...
                RU_valid,InspBU_valid,ObsYears_valid,Ncurve,OptBoundsData,...
                SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);
            TrainSmootherRun();          
        catch
            disp(' ')
            disp('Constraints Failiur')
            LogLikValidation=-99^10;
            loglik=-99^10;
            break;
        end
        % Stopping Criteria
        SumLL(RegressLoop,1)=gather(sum(LogLikValidation(:)));
        if RegressLoop>1
            if abs(SumLL(RegressLoop-1,1))-abs(SumLL(RegressLoop,1))<TolRegress
                Converge=1;
                RegressLoop=MultiPass;
                InirilizedEx=InirilizedExStore;
                InirilizedVar=InirilizedVarStore;
                LogLikValidation=LoglikStoreValidation;
            else
                InirilizedExStore=InirilizedEx;
                InirilizedVarStore=InirilizedVar;
                LoglikStoreValidation=LogLikValidation;
            end
        else
            InirilizedExStore=InirilizedEx;
            InirilizedVarStore=InirilizedVar;
            LoglikStoreValidation=LogLikValidation;
        end
        RegressLoop=RegressLoop+1;
    end
    %% Evaluate the estimated parameters
    if RegressLoop>=MultiPass || Converge || args{6}==2
    % initial speed state
    Converge=1;
    init_x(2,:)=AKr*InirilizedEx;
    init_V(2,2,:)=diag(AKr*InirilizedVar*AKr') + (Sigma_W0^2*ones(size(AKr,1),1));%(param(4).^2).*(DifferenceObs)+(param(6).^2);

    % Evaluate the Loglikelihood
    [x, Var, ~, loglik,~,~] = kalman_filter(y, A, C, Q, R, Re...
        , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,...
        ObsYears,Ncurve,OptBoundsData,SpeedConstraints(1,1),...
        SpeedConstraints,GPUCompute,Nsigma);
    TrainSmootherRun(); 
    end
    
    if ~isempty(TrueSpeedValues)
        figure(1)
        plot(TrueSpeedValues,init_x(2,:),'o')
        axis equal
        hold on
        plot([0,-.6],[0,-.6]);
        xlim([-0.6 0])
        ylim([-0.6 0])
        ylabel('True Speed')
        xlabel('Estiamted Speed')
        grid on
        hold off
    end
end


% %     Constrained Kernel Regression
%     Nval=length(InirilizedEx);
%     AKr=Nval.*AKr;
%     Amat=-2.*AKr.*InirilizedEx(:,1)';
%     bmat=1.5*ones(length(Amat),1);
%     Aeq=ones(1,Nval);
%     beq=1;
%     Hmat=eye(Nval);
%     fmat=zeros(Nval,1);
%     Pu = quadprog(Hmat,fmat,Amat,bmat,Aeq,beq);
%     AKr = AKr.*Pu';

%     init_V(2,2,:)=diag(AKr*InirilizedVar(:,:,1)*AKr' + diag(Sigma_W0^2*ones(size(AKr,1),1)));