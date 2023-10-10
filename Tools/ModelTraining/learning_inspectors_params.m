inspectors_params_updated = 0;
while (LLcr/LLprev)<= StopCr && (OI  || ( StallInit1<StallVal1))
    LLprev=LLcr;
    for model_i=1:length(app.Tree.Children) % model type iterations: 2: SSM-KR or 1: SSM
        for elem_j=1:length(app.Tree.Children(model_i).Children)
            % Read data
            ElementName= app.Tree.Children(model_i).Children(elem_j).Text;
            fprintf('\nEstimating inspector parameters of %s',ElementName);
            app.Tree.SelectedNodes=app.Tree.Children(model_i).Children(elem_j);
            ElementData=load([FullPath 'TrainingData_' erase(ElementName,"/") '_' sprintf('%d',Ncurve) '.mat']);
            ElementData=struct2cell(ElementData);
            ElementData=ElementData{1};
            Index=find(strcmp(app.AllElementsParameters(:,1),ElementName));
            % use the estimated inspectors parameters as initial guess for the
            % next set of elements
            if LLcr==-10^12
                stored_inspectors_params = [stored_inspectors_params; AllElementsParameters{Index,3}];
                stored_init_estim = [stored_init_estim; zeros(size(stored_inspectors_params,1),4)];
            end
            if inspectors_params_updated
                for ej=1:length(InspectorsData{1}(:,1))
                    ind_elem_j = find(stored_inspectors_params(:,1) == InspectorsData{1}(ej,1));
                    if ~isempty(ind_elem_j)
                        stored_inspectors_params(min(ind_elem_j),2:3) = InspectorsData{1}(ej,2:3);
                        stored_init_estim(min(ind_elem_j),1) = 0;
                        stored_init_estim(min(ind_elem_j),2) = sqrt(cov_params(4,4,ej));
                        stored_init_estim(min(ind_elem_j),3) = PARAM(2);
                        stored_init_estim(min(ind_elem_j),4) = sqrt(cov_params(5,5,ej));
                    end
                end
            end
            [~, iau] = unique(stored_inspectors_params(:,1));
            stored_inspectors_params = stored_inspectors_params(iau,:);
            stored_init_estim = stored_init_estim(iau,:);
            if inspectors_params_updated
                [~, ~,iai] = intersect(AllElementsParameters{Index,3}(:,1),stored_inspectors_params(:,1));
                AllElementsParameters{Index,3} = stored_inspectors_params(iai,:);
                init_estim = stored_init_estim(iai,:);
                check_init = find(init_estim(:,2)==0);
                if ~isempty(check_init)
                    init_estim(check_init,2) = 1;
                    init_estim(check_init,3) = PARAM(2);
                    init_estim(check_init,4) = 12;
                end
            else
                %init_estim = [ 0, 1, PARAM(2), 12];
                [~, ~,iai] = intersect(AllElementsParameters{Index,3}(:,1),stored_inspectors_params(:,1));
                init_estim = stored_init_estim(iai,:);
                init_estim(:,1) = 0;
                init_estim(:,2) = 1;
                init_estim(:,3) = PARAM(2);
                init_estim(:,4) = 12;
            end
            InspectorsData=load([FullPath 'Inspectors_' erase(ElementName,"/") '.mat']);
            InspectorsData=struct2cell(InspectorsData);
            InspectorsData{1} = AllElementsParameters{Index,3};
            PARAM = AllElementsParameters{Index,2};
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
            inspectors_params_updated = 1;
            %InspectorsData{1} = AllElementsParameters{Index,3};
        end
    end
    for ej=1:length(InspectorsData{1}(:,1))
        ind_elem_j = find(stored_inspectors_params(:,1) == InspectorsData{1}(ej,1));
        if ~isempty(ind_elem_j)
            stored_inspectors_params(min(ind_elem_j),2:3) = InspectorsData{1}(ej,2:3);
        end
    end
    
    for model_i=1:length(app.Tree.Children) % model type iterations: 2: SSM-KR or 1: SSM
        for elem_j=1:length(app.Tree.Children(model_i).Children)
            % Read data
            ElementName= app.Tree.Children(model_i).Children(elem_j).Text;
            app.Tree.SelectedNodes=app.Tree.Children(model_i).Children(elem_j);
            ElementData=load([FullPath 'TrainingData_' erase(ElementName,"/") '_' sprintf('%d',Ncurve) '.mat']);
            ElementData=struct2cell(ElementData);
            ElementData=ElementData{1};
            Index=find(strcmp(app.AllElementsParameters(:,1),ElementName));
            [~, ~,iai] = intersect(AllElementsParameters{Index,3}(:,1),stored_inspectors_params(:,1));
            AllElementsParameters{Index,3} = stored_inspectors_params(iai,:);
            InspectorsData{1} = AllElementsParameters{Index,3};
            PARAM = AllElementsParameters{Index,2};
            if model_i==1
                IncludeStructuralAtt=0;
            else
                IncludeStructuralAtt=1;
                NumAttributes = app.AllElementsParameters{Index,2}{5} ;
            end
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
            LLcr_cat = LLcr_cat + sum(LogLikVal);
        end
    end
    LLcr = LLcr_cat;
    LLcr_cat = 0;
    inspectors_params_updated = 0;
    if LLcr<LLprev
        stored_inspectors_params = Stored_InspectorData;
        break
    else
        Stored_InspectorData = stored_inspectors_params;
    end
    if LLcr/LLprev>0.95
        StallInit1=StallInit1+1;
    end
end
StallInit1=0;
