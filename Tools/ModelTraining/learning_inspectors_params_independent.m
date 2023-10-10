inspectors_params_updated = 0 * ones(1,total_num_elements);

while (LLcr_cats(j)/LLprev_cats(j))<= StopCr && (OI  || ( StallInit1<StallVal1))
    LLprev_cats(j) = LLcr_cats(j);
    % use the estimated inspectors parameters as initial guess for the
    % next set of elements
    
    PARAM = AllElementsParameters{Index,2};
    if LLcr_cats(j) == -10^12
        stored_inspectors_params = AllElementsParameters{Index,3};
        stored_init_estim = zeros(size(stored_inspectors_params,1),4);
    end
    init_estim = stored_init_estim;
    init_estim(:,1) = 0;
    init_estim(:,2) = 1;
    init_estim(:,3) = PARAM(2);
    init_estim(:,4) = 12;
    InspectorsData=load([FullPath 'Inspectors_' erase(ElementName,"/") '.mat']);
    InspectorsData=struct2cell(InspectorsData);
    InspectorsData{1} = AllElementsParameters{Index,3};
    
    if model_i==1
        IncludeStructuralAtt=0;
    else
        IncludeStructuralAtt=1;
        NumAttributes = app.AllElementsParameters{Index,2}{5} ;
    end
    if IncludeStructuralAtt
        AnnModel = AllElementsParameters{Index,4};
        AnnModel.categ_st_att_ind = Material_index(elem_j);
    end
    
    % pre-processing data
    Struc_Att_preprocessing();
    
    % learning inspectors parameters
    OnlineInference_Main();
    
    AllElementsParameters{Index,3} = InspectorsData{1};
    inspectors_params_updated(j) = 1;
    %InspectorsData{1} = AllElementsParameters{Index,3};
    % Read data
    if IncludeStructuralAtt
        AnnModel = AllElementsParameters{Index,4};
        AnnModel.categ_st_att_ind = Material_index(elem_j);
        % validation with validation set
        if ~isempty(AnnModel)
            if isempty(AnnModel.theta_nn)
                StrucAttributes=0;
                StrucAttributes_Test=0;
            else
                StrucAttributes=ElementData.ModelValid.StrucAtt;
                StrucAttributes_Test=ElementData.ModelTest.StrucAtt;
            end
        else
            StrucAttributes=ElementData.ModelValid.StrucAtt;
            StrucAttributes_Test=ElementData.ModelTest.StrucAtt;
        end
        [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},...
            [],Ncurve,OptimizationProceedure,OptLevel,SpeedConstraints,...
            1,StrucAttributes,KRparam(2:end),...
            KRparam(1),KernelParameters{1},KernelParameters{2},...
            InitialEx,InitialVar,[],AnnModel);
        [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            StrucAttributes_Test,KRparam(2:end),KRparam(1),...
            KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[],AnnModel);
        TestLogLik=sum(LogLikVal_Test);
        fprintf('Valid L.L. : %d \n',sum(LogLikVal))
        fprintf('Test L.L. : %d \n',TestLogLik)
    else
        % validation with validation set
        [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1);
        [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,...
            SpeedConstraints,1);
        TestLogLik=sum(LogLikVal_Test);
        fprintf('Valid L.L. : %d \n',sum(LogLikVal))
        fprintf('Test L.L. : %d \n',TestLogLik)
    end
    plot(app.UIAxes,datetime('now'),gather(sum(TestLogLik)),'o');
    hold(app.UIAxes,'on');
    grid(app.UIAxes,'on');
    box(app.UIAxes,'on');
    ylabel(app.UIAxes,'Test set Logliklihood')
    xlabel(app.UIAxes,'Time')
    drawnow;
    LLcr_cats(j) = sum(LogLikVal);
    if LLcr_cats(j) < LLprev_cats(j)
        AllElementsParameters{Index,3} = stored_inspectors_params;
        InspectorsData{1} = stored_inspectors_params;
        inspectors_params_updated(j) = 0;
        break
    else
        stored_inspectors_params = AllElementsParameters{Index,3};
    end
    if LLcr_cats(j) / LLprev_cats(j)>0.95
        StallInit1=StallInit1+1;
    end
end
StallInit1=0;
