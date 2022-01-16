OptLevel=2;
StopCr=10^-3;
LLcr=-10^12;
LLprev=-10^16;
j=1;

StallVal1=10;
StallVal2=10;
StallInit1=0;
StallInit2=0;


while (LLcr-LLprev)>StopCr && StallInit2<StallVal2
    while (LLcr-LLprev)>StopCr && StallInit1<StallVal1
        LLprev=LLcr;
        for i=1:1%length(InspectorsData{1}(:,1))
            sprintf('Inspector Num: %d /%d',i,length(InspectorsData{1}(:,1)))
            InspectParam=[UpdatedInspectorsData{1}(i,end) UpdatedInspectorsData{1}(i,2)];
            InspectID=InspectorsData{1}(i,1);
            if IncludeStructuralAtt
                [Engparam,~,~,fx_NR]=Newton_Raphson_par_one(@(InspectParam) AnalysisObjective(ElementData,...
                    InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                    Q, x0, s2_X0,PARAM,InspectParam,InspectorsData{1},...
                    InspectID,Ncurve,OptimizationProceedure,OptLevel,SpeedConstraints,...
                    1,StructuralAttributes,KRparam(2:end),KRparam(1),KernelParameters{1},...
                    KernelParameters{2},InitialEx,InitialVar,[]),InspectParam,'log_transform','no',...
                    'output','original','laplace','no','convergence_tol',1E-4,...
                    'bounds',bounds,'nb_failed_iteration_limit',1);
            else
                [Engparam,~,~,fx_NR]=Newton_Raphson_par_one(@(InspectParam) AnalysisObjective(ElementData,...
                    InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                    Q, x0, s2_X0,PARAM, InspectParam, InspectorsData{1},...
                    InspectID,Ncurve,OptimizationProceedure,OptLevel,SpeedConstraints,...
                    1),InspectParam,'log_transform','no',...
                    'output','original','laplace','no','convergence_tol',1E-4,...
                    'bounds',bounds,'nb_failed_iteration_limit',1);
            end
            InspectorsData{1}(i,3) = Engparam(1);
            InspectorsData{1}(i,2) = Engparam(2);
        end
        
        if IncludeStructuralAtt
            % validation with validation set
            if isempty(InitialEx)
                StrucAttributes=0;
                StrucAttributes_Test=0;
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
                InitialEx,InitialVar,[]);
            [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
                InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                OptimizationProceedure,OptLevel,SpeedConstraints,1,...
                StrucAttributes_Test,KRparam(2:end),KRparam(1),...
                KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
            TestLogLik=sum(LogLikVal_Test)
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
            TestLogLik=sum(LogLikVal_Test)
        end
        plot(app.UIAxes,datetime('now'),gather(sum(TestLogLik)),'o');
        hold(app.UIAxes,'on');
        grid(app.UIAxes,'on');
        box(app.UIAxes,'on');
        ylabel(app.UIAxes,'Test set Logliklihood')
        xlabel(app.UIAxes,'Time')
        pause(0.01)
        LLcr=sum(LogLikVal);
        if LLcr<LLprev
            break
        end
        if LLcr/LLprev>0.95
            StallInit1=StallInit1+1;
        end
    end
    StallInit1=0;
    
    OptLevel=3;
    j=j+1;
    if IncludeStructuralAtt
        [Qparam,~,~,~]=Newton_Raphson_par(@(PARAM) AnalysisObjective(ElementData,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            StructuralAttributes,KRparam(2:end),KRparam(1),...
            KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]),...
            PARAM,'log_transform','no','output','original','laplace',...
            'no','convergence_tol',1E-4,'bounds',boundsQ);
        PARAM=Qparam;
        StructuralAttributes=ElementData.StrucAtt;
        [KRparam,~,~,~]=Newton_Raphson_par(@(Kr_Param) AnalysisObjective(ElementData,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
            x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            ElementData.StrucAtt,Kr_Param(2:end),Kr_Param(1),...
            KernelParameters{1},KernelParameters{2},[],[],[])...
            ,Kr_Param,'log_transform','no','output','original','laplace','no',...
            'convergence_tol',1E-4,'bounds',boundKr,'delta',1E-3);
        Kr_Param=KRparam;
        InirilizedEx=ones(KernelParameters{2},1)*0;
        [fx_NR,~,InitialEx,InitialVar,X_ControlPoints]=AnalysisObjective(ElementData,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
            x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            ElementData.StrucAtt,KRparam(2:end),KRparam(1),KernelParameters{1},...
            KernelParameters{2},InirilizedEx,[],[]);
        RegressionModel.InirilizedEx=InitialEx;
        RegressionModel.InirilizedVar=InitialVar;
        RegressionModel.Sigma_W0=KRparam(1);
        RegressionModel.KernelType=KernelParameters{1};
        RegressionModel.X_ControlPoints=X_ControlPoints;
        RegressionModel.Kernel_l=KRparam(2:end);
        
        % validation with validation set
        [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            ElementData.ModelValid.StrucAtt,KRparam(2:end),KRparam(1),...
            KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
        % Testingwith Testing set
        [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
            ElementData.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
            KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
        TestLogLik=sum(LogLikVal_Test)
        
    else
        [Qparam,~,~,fx_NR]=Newton_Raphson_par(@(PARAM) AnalysisObjective(ElementData,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,...
            SpeedConstraints,1),PARAM,'log_transform','no','output','original','laplace',...
            'no','convergence_tol',1E-4,'bounds',boundsQ);
        PARAM=Qparam;
        % validation with validation set
        [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1);
        [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1);
        TestLogLik=sum(LogLikVal_Test)
    end
    
    hold(app.UIAxes,'on');
    plot(app.UIAxes,datetime('now'),gather(TestLogLik),'o');
    grid(app.UIAxes,'on');
    box(app.UIAxes,'on');
    ylabel(app.UIAxes,'Test Set Logliklihood')
    xlabel(app.UIAxes,'Time')
    pause(0.1);
    
    % save parameters
    
    AllElementsParameters=app.AllElementsParameters(:,1);
    AllElementsParameters{Index,2}=PARAM;
    AllElementsParameters{Index,3}=InspectorsData{1};
    AllElementsParameters{Index,4}=RegressionModel;
    ElementMetaData=load([FullPath 'MetaData_' ElementName  '.mat']);
    ElementMetaData=struct2cell(ElementMetaData);
    ElementMetaData=ElementMetaData{1};
    if i==1
        AttOrder=load([FullPath 'RegressAtt_' ElementName  '.mat']);
        AttOrder=struct2cell(AttOrder);
        AttOrder=AttOrder{1};
        AttOrder=str2double(cellstr(AttOrder));
        AttOrder=AttOrder(1:NumAttributes)+4;
    else
        AttOrder=[];
    end
    AllElementsParameters{Index,5}=AttOrder;
    AllElementsParameters{Index,9}=ElementMetaData;
    AllElementsParameters{Index,10}=Ncurve;
    save([SavePath '/AutoSave_' ElementName '.mat'],'AllElementsParameters');
    
    
    LLcr=sum(LogLikVal);
    if LLcr<LLprev
        break
    end
    if LLcr/LLprev>0.95
        StallInit2=StallInit2+1;
    end
    
end
%% Testing with Testing set
if IncludeStructuralAtt
    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelTest,...
        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
        OptimizationProceedure,OptLevel,...
        SpeedConstraints,1,ElementData.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
        KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
    TestLogLik=sum(LogLikVal)
else
    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelTest,...
        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
        OptimizationProceedure,OptLevel,...
        SpeedConstraints,1);
    TestLogLik=sum(LogLikVal)
end





