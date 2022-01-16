% prior assumptions about an inspector. 0: no prior
PriorKnowledgeChk=0;MeanPrior=0;StdPrior=0;
OptLevel=2;
StopCr2= 0.999;   %10^-3;
StopCr1 = 0.999;
LLcr=-10^12;
LLprev=-10^18;
j=1;

if get(app.LogParams,'Value')==1 && isempty(ModelParameters)
    save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/UpdatedInspectorsData%s%s.mat',pwd,Data_filename,date),'UpdatedInspectorsData');
end

StallVal1=10;
StallVal2=10;
StallInit1=0;
StallInit2=0;



firstrun = 0;

if OptimizationAlgorithmIndex==1
    OI = 0;
else
    OI = 1;
end   
epo = 1;
indice = 0;    
while  (LLcr/LLprev)<= StopCr2 &&  (OI  || ( StallInit2<StallVal2))
    while (LLcr/LLprev)<= StopCr1 &&  (OI  || ( StallInit1<StallVal1))
  
        LLprev=LLcr;
        if OptimizationAlgorithmIndex==1
            for i=1:length(UpdatedInspectorsData{1}(:,1))
                fprintf('Inspector Num: %d /%d',i,length(UpdatedInspectorsData{1}(:,1)))
                if OperationIndex>=3
                    InspectParam=[UpdatedInspectorsData{1}(i,end) UpdatedInspectorsData{1}(i,2)];
                else
                    InspectParam= UpdatedInspectorsData{1}(i,end);
                end
                    
                InspectID=UpdatedInspectorsData{1}(i,1);
                if IncludeStructuralAtt
                    [Engparam,~,~,fx_NR]=Newton_Raphson_par_one(@(InspectParam) AnalysisObjective(MdataEngy,...
                        UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,InspectParam,UpdatedInspectorsData{1},...
                        InspectID,NTr,DataTransformedSpace,OptimizationProceedure,...
                        ParComp,OptLevel,GlobalCondData,PriorKnowledgeChk, MeanPrior,...
                        StdPrior,1,GradCompute,StructuralAttributes,KRparam(2:end),KRparam(1),KernelParameters{1},...
                        KernelParameters{2},InitialEx,InitialVar,[]),InspectParam,'log_transform','no',...
                        'output','original','laplace','no','convergence_tol',1E-4,...
                        'bounds',bounds,'nb_failed_iteration_limit',1);
                else
                    [Engparam,~,~,fx_NR]=Newton_Raphson_par_one(@(InspectParam) AnalysisObjective(MdataEngy,...
                        UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,InspectParam,UpdatedInspectorsData{1},...
                        InspectID,NTr,DataTransformedSpace,OptimizationProceedure,...
                        ParComp,OptLevel,GlobalCondData,PriorKnowledgeChk, MeanPrior,...
                        StdPrior,1,GradCompute),InspectParam,'log_transform','no',...
                        'output','original','laplace','no','convergence_tol',1E-4,...
                        'bounds',bounds,'nb_failed_iteration_limit',1);
                end
                UpdatedInspectorsData{1}(i,3)=Engparam(1); % inspector std
                if OperationIndex>=3
                    UpdatedInspectorsData{1}(i,2) = Engparam(2); % inspector bias
                else
                    UpdatedInspectorsData{1}(i,2) = 0;           % inspector bias
                end
            end
            
        else
            OnlineInference_Main()
            StallInit2 = 1 + StallVal2;
            epo = 1;
        end
        
        if IncludeStructuralAtt
            % validation with validation set
            if isempty(InitialEx)
                StrucAttributes=0;
                StrucAttributes_Test=0;
            else
                StrucAttributes=MdataEngy.ModelValid.StrucAtt;
                StrucAttributes_Test=MdataEngy.ModelTest.StrucAtt;
            end
            [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelValid,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},...
                [],NTr,DataTransformedSpace,OptimizationProceedure,...
                ParComp,OptLevel,GlobalCondData,PriorKnowledgeChk, MeanPrior,...
                StdPrior,1,GradCompute,StrucAttributes,KRparam(2:end),...
                KRparam(1),KernelParameters{1},KernelParameters{2},...
                InitialEx,InitialVar,[]);
            [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(MdataEngy.ModelTest,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute,StrucAttributes_Test,KRparam(2:end),KRparam(1),...
                KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
            TestLogLik=sum(LogLikVal_Test);
            fprintf('Test L.L. : %d \n',TestLogLik)
        else
            % validation with validation set
            [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelValid,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute);
            [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(MdataEngy.ModelTest,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute);
            TestLogLik=sum(LogLikVal_Test);
            fprintf('Test L.L. : %d \n',TestLogLik)
        end
        % Check the inspector parameters
        if app.DTransformedSpace.Value && ~isempty(app.TrueInspectorsData) && app.LogParams_2.Value
            cla(app.OutputPlot);
            xlim(app.OutputPlot,[1 7]);
            ylim(app.OutputPlot,[1 7]);
            plot(app.OutputPlot,[1,7],[1,7]);
            hold(app.OutputPlot,'on');
            plot(app.OutputPlot,UpdatedInspectorsData{1}(:,end),app.TrueInspectorsData{1}(:,1),'o');
            hold(app.OutputPlot,'off');
            grid(app.OutputPlot,'on')
            box(app.OutputPlot,'on');
            ylabel(app.OutputPlot,'Inspector True $\sigma(I_i)$','interpreter','latex')
            xlabel(app.OutputPlot,'Inspector Estimated $\sigma(I_i)$','interpreter','latex')
            drawnow;
        elseif app.LogParams_2.Value
            histogram(app.OutputPlot,UpdatedInspectorsData{1}(:,end));
            xlabel(app.OutputPlot,'Histogram of Inspector Estimated Params $\sigma(I_i)$','interpreter','latex');
            drawnow;
        end
        if get(app.LoglikelihoodPlot,'value')==1
            plot(app.OutputPlot1,datetime('now'),gather(sum(TestLogLik)),'o');
            hold(app.OutputPlot1,'on');
            grid(app.OutputPlot1,'on');
            box(app.OutputPlot1,'on');
            ylabel(app.OutputPlot1,'Test set Logliklihood')
            xlabel(app.OutputPlot1,'Time')
            drawnow;
        end
        LLcr=sum(LogLikVal)
        if LLcr<LLprev
            UpdatedInspectorsData = Stored_UpdatedInspectorData;
            break
        else 
            Stored_UpdatedInspectorData = UpdatedInspectorsData;
        end
        if get(app.LogParams,'Value')==1
            save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/UpdatedInspectorsData%s%s.mat',pwd,Data_filename,date),'UpdatedInspectorsData');
        end
        if LLcr/LLprev>0.95
            StallInit1=StallInit1+1;
        end
    end
    StallInit1=0;
    if isempty(find(OptimizationProceedure==[3,5,7]))
        OptLevel=3;
        j=j+1;
        if IncludeStructuralAtt
            [Qparam,~,~,~]=Newton_Raphson_par(@(PARAM) AnalysisObjective(MdataEngy,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute,StructuralAttributes,KRparam(2:end),KRparam(1),...
                KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]),...
                PARAM,'log_transform','no','output','original','laplace',...
                'no','convergence_tol',1E-4,'bounds',boundsQ);
            PARAM=Qparam;
            StructuralAttributes=MdataEngy.StrucAtt;
            [KRparam,~,~,~]=Newton_Raphson_par(@(Kr_Param) AnalysisObjective(MdataEngy,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
                x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,GradCompute,...
                MdataEngy.StrucAtt,Kr_Param(2:end),Kr_Param(1),...
                KernelParameters{1},KernelParameters{2},[],[],[])...
                ,Kr_Param,'log_transform','no','output','original','laplace','no',...
                'convergence_tol',1E-4,'bounds',boundKr,'delta',1E-3);
            Kr_Param=KRparam;
            InirilizedEx=ones(KernelParameters{2},1)*0;
            [fx_NR,~,InitialEx,InitialVar,X_ControlPoints]=AnalysisObjective(MdataEngy,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
                x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,GradCompute,...
                MdataEngy.StrucAtt,KRparam(2:end),KRparam(1),KernelParameters{1},...
                KernelParameters{2},InirilizedEx,[],TrueSpeedValues);
            RegressionModel.InirilizedEx=InitialEx;
            RegressionModel.InirilizedVar=InitialVar;
            RegressionModel.Sigma_W0=KRparam(1);
            RegressionModel.KernelType=KernelParameters{1};
            RegressionModel.X_ControlPoints=X_ControlPoints;
            RegressionModel.Kernel_l=KRparam(2:end);
            if get(app.LogParams,'value')==1
                save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/KR_PARAM%s%s.mat',pwd,Data_filename,date),'RegressionModel');
            end
            % validation with validation set
            [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelValid,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute,MdataEngy.ModelValid.StrucAtt,KRparam(2:end),KRparam(1),...
                KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
            % Testingwith Testing set
            [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(MdataEngy.ModelTest,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute,MdataEngy.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
                KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
            TestLogLik=sum(LogLikVal_Test);
            fprintf('Test L.L. : %d \n',TestLogLik)
        else
            [Qparam,~,~,fx_NR]=Newton_Raphson_par(@(PARAM) AnalysisObjective(MdataEngy,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute),PARAM,'log_transform','no','output','original','laplace',...
                'no','convergence_tol',1E-4,'bounds',boundsQ);
            PARAM=Qparam;
            % validation with validation set
            [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelValid,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute);
            [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(MdataEngy.ModelTest,...
                UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
                DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
                GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
                GradCompute);
            TestLogLik=sum(LogLikVal_Test);
            fprintf('Test L.L. : %d \n',TestLogLik)
        end
        if get(app.LoglikelihoodPlot,'value')==1
            hold(app.OutputPlot1,'on');
            plot(app.OutputPlot1,datetime('now'),gather(TestLogLik),'o');
            grid(app.OutputPlot1,'on');
            box(app.OutputPlot1,'on');
            ylabel(app.OutputPlot1,'Test Set Logliklihood')
            xlabel(app.OutputPlot1,'Time')
            drawnow;
        end
    end
    LLcr=sum(LogLikVal);
    if LLcr<LLprev
        RegressionModel = Stored_RegressionModel;
        Qparam = Stored_QParam;
        PARAM = Qparam;
        break
    else 
        Stored_RegressionModel = RegressionModel;
        Stored_QParam = Qparam;
    end
    if get(app.LogParams,'Value')==1
        save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/UpdatedInspectorsData%s%s.mat',pwd,Data_filename,date),'UpdatedInspectorsData');
        save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/PARAM%s%s.mat',pwd,Data_filename,date),'Qparam');
        if IncludeStructuralAtt
            save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/KR_PARAM%s%s.mat',pwd,Data_filename,date),'RegressionModel');
        end
    end
    if LLcr/LLprev>0.95
        StallInit2=StallInit2+1;
    end
end
%% Testing with Testing set
if IncludeStructuralAtt
    [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelTest,...
        UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
        Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
        DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
        GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
        GradCompute,MdataEngy.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
        KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]);
    TestLogLik=sum(LogLikVal);
    fprintf('Test L.L. : %d\n',TestLogLik)
    save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/KR_PARAM%s%s.mat',pwd,Data_filename,date),'RegressionModel');
else
    [~,~,~,~,~,LogLikVal]=AnalysisObjective(MdataEngy.ModelTest,...
        UpdatedInspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
        Q, x0, s2_X0,PARAM,[0 0],UpdatedInspectorsData{1},[],NTr,...
        DataTransformedSpace,OptimizationProceedure,ParComp,OptLevel,...
        GlobalCondData,PriorKnowledgeChk, MeanPrior, StdPrior,1,...
        GradCompute);
    TestLogLik=sum(LogLikVal);
    fprintf('Test L.L. : %d\n',TestLogLik)
end
save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/UpdatedInspectorsData%s%s.mat',pwd,Data_filename,date),'UpdatedInspectorsData');
save(sprintf('%s/Tools/ModelTraining_Generic/ParametersLog/PARAM%s%s.mat',pwd,Data_filename,date),'Qparam');
LogLik=LLcr;  