OptLevel=2;
total_num_elements =(length(app.Tree.Children(1).Children)+length(app.Tree.Children(2).Children));
StopCr=0.999;
LLcr=-10^12;
LLcr_cat = 0;
LLprev=-10^18;
StallVal1=3;
StallInit1=0;


LLcr_cats = -10^12 * ones(1,total_num_elements);
LLprev_cats = -10^18 * ones(1,total_num_elements);
StallVal2 = 3 * ones(1,total_num_elements);
StallInit2 = 0 * ones(1,total_num_elements);

Stored_RegressionModel = cell(length(app.Tree.Children),max(length(app.Tree.Children(1).Children),length(app.Tree.Children(2).Children)));
Stored_QParam = cell(length(app.Tree.Children),max(length(app.Tree.Children(1).Children),length(app.Tree.Children(2).Children)));
stored_inspectors_params = [];
stored_init_estim = [];

keep_learning = 1;

OI = 1;
bnn_available = zeros(1,total_num_elements);

while keep_learning
    
       if joint_inspector_train
            learning_inspectors_params(); %-oops
       end
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    OptLevel=3;
    j=1;
    for model_i=1:length(app.Tree.Children) % model type iterations: 2: SSM-KR or 1: SSM
        for elem_j=1:length(app.Tree.Children(model_i).Children)
            loop_over_category = 1; % to differentiate independent from joint estimation 
            while loop_over_category
                ElementName= app.Tree.Children(model_i).Children(elem_j).Text;
                d.Value = (elem_j/NumCats)/1.2;
                d.Message=sprintf('\nStarted training structural category %s',ElementName);
                disp("====================================================================================")
                disp(d.Message)
                if d.CancelRequested
                    break
                end
                pause(0.001);
                % Read data
                app.Tree.SelectedNodes=app.Tree.Children(model_i).Children(elem_j);
                ElementData=load([FullPath 'TrainingData_' erase(ElementName,"/") '_' sprintf('%d',Ncurve) '.mat']);
                ElementData=struct2cell(ElementData);
                ElementData=ElementData{1};     
                Index=find(strcmp(app.AllElementsParameters(:,1),ElementName));
                if model_i==1
                    IncludeStructuralAtt=0;
                else
                    IncludeStructuralAtt=1;
                    NumAttributes = app.AllElementsParameters{Index,2}{5} ;
                end
                Struc_Att_preprocessing();
                if joint_inspector_train
                    [~, ~, iai] = intersect(AllElementsParameters{Index,3}(:,1),stored_inspectors_params(:,1));%-oops
                    AllElementsParameters{Index,3} = stored_inspectors_params(iai,:);%-oops
                end
                InspectorsData{1} = AllElementsParameters{Index,3};
                PARAM = AllElementsParameters{Index,2};
                if ~joint_inspector_train
                    learning_inspectors_params_independent(); %-oops
                end
                if IncludeStructuralAtt
                    AnnModel = AllElementsParameters{Index,4};
                    AnnModel.categ_st_att_ind = Material_index(elem_j);
                    if ~stopped_learning_1(elem_j)
                        [Qparam,~,~,~]=Newton_Raphson(@(PARAM) AnalysisObjective(ElementData,...
                            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
                            StructuralAttributes,KRparam(2:end),KRparam(1),...
                            KernelParameters{1},KernelParameters{2},[],[],[],AnnModel),...
                            PARAM,'log_transform','no','output','original','laplace',...
                            'no','convergence_tol',5E-1,'bounds',boundsQ);
                        PARAM = Qparam;
                        StructuralAttributes=ElementData.StrucAtt;
                        AnnModel.evaluate_mode = 0;
                        [fx_NR,~,~,~,~,~, AnnModel_params]=AnalysisObjective(ElementData,...
                            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
                            x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                            OptimizationProceedure,OptLevel,SpeedConstraints,1,...
                            ElementData.StrucAtt,KRparam(2:end),KRparam(1),...
                            KernelParameters{1},KernelParameters{2},[],[],[],AnnModel);
                        AnnModel = AnnModel_params;
                        bnn_available(elem_j) = 1;
                    end
                    AnnModel.cont_training = 1;
                    AnnModel.evaluate_mode = 1; % ensures AnalysisObjective returns loglik using TAGI to initialize speeds
                    
                    % validation with validation set
                    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,SpeedConstraints,1,...
                        ElementData.ModelValid.StrucAtt,KRparam(2:end),KRparam(1),...
                        KernelParameters{1},KernelParameters{2},[],[],[],AnnModel);
                    % Testingwith Testing set
                    [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,SpeedConstraints,1,...
                        ElementData.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
                        KernelParameters{1},KernelParameters{2},[],[],[],AnnModel);
                    TestLogLik=sum(LogLikVal_Test);
                    ValLogLik=sum(LogLikVal);
                    fprintf('Valid L.L. : %d \n',ValLogLik)
                    fprintf('Test L.L. : %d \n',TestLogLik)
                else
                    if ~stopped_learning_1(elem_j)
                        %% commented this temp - note that convergence tol is relaxed from 1e-4 to 1e-1
                        [Qparam,~,~,fx_NR]=Newton_Raphson(@(PARAM) AnalysisObjective(ElementData,...
                            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                            Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                            OptimizationProceedure,OptLevel,...
                            SpeedConstraints,1),PARAM,'log_transform','no','output','original','laplace',...
                            'no','convergence_tol',5E-1,'bounds',boundsQ);
                        PARAM = Qparam;
%                         Qparam = PARAM;
                end
%                     Qparam=PARAM;
                    % validation with validation set
                    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelValid,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,SpeedConstraints,1);
                    [~,~,~,~,~,LogLikVal_Test]=AnalysisObjective(ElementData.ModelTest,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,SpeedConstraints,1);
                    TestLogLik=sum(LogLikVal_Test);
                    ValLogLik=sum(LogLikVal);
                    fprintf('Valid L.L. : %d \n',ValLogLik)
                    fprintf('Test L.L. : %d \n',TestLogLik)
                end
                
                hold(app.UIAxes,'on');
                plot(app.UIAxes,datetime('now'),gather(TestLogLik),'o');
                grid(app.UIAxes,'on');
                box(app.UIAxes,'on');
                ylabel(app.UIAxes,'Test Set Logliklihood')
                xlabel(app.UIAxes,'Time')
                drawnow;
                
                LLcr_cats(j) = sum(LogLikVal);
                
                % save parameters
                
                AllElementsParameters{Index,2}=PARAM;
                AllElementsParameters{Index,3}=InspectorsData{1};
                AllElementsParameters{Index,4}=AnnModel;
                if strcmp(app.Switch.Value, 'Real')
                    ElementMetaData=load([FullPath 'MetaData_' erase(ElementName,"/")  '.mat']);
                    ElementMetaData=struct2cell(ElementMetaData);
                    ElementMetaData=ElementMetaData{1};
                    if IncludeStructuralAtt
                        % reasoning is @ RunTraining_unit() because some
                        % structural attributes are constant for elements
                        % therefore they are reomved, in this instance it is
                        % the material
                        if sum(strcmp(ElementName,{'Bande médiane', 'Tirants', 'Élément en élastomère', 'Toiture'})) == 1
                            AttOrder=(1:NumAttributes)+5;
                        else 
                            AttOrder=(1:NumAttributes)+4;
                        end
                    else
                        AttOrder=[];
                    end
                    AllElementsParameters{Index,5}=AttOrder;
                    AllElementsParameters{Index,9}=ElementMetaData;
                end
                AllElementsParameters{Index,10}=Ncurve;
                
                if LLcr_cats(j)<LLprev_cats(j)
                    if IncludeStructuralAtt
                        if ~isempty(Stored_RegressionModel{model_i, elem_j})
                            AnnModel = Stored_RegressionModel{model_i, elem_j};
                        end
                    end
                    if ~isempty(Stored_QParam{model_i, elem_j})
                        Qparam = Stored_QParam{model_i, elem_j};
                    end
                    PARAM = Qparam;
                    AllElementsParameters{Index,2}=PARAM;
                    AllElementsParameters{Index,4}=AnnModel;
                    stopped_learning_1(elem_j) = 1;
                    save([SavePath '/AutoSave_' 'learning_status_1.mat'],'stopped_learning_1');
                else
                    if IncludeStructuralAtt
                        Stored_RegressionModel{model_i, elem_j} = AnnModel;
                    end
                    Stored_QParam{model_i, elem_j} = AllElementsParameters{Index,2};%Qparam;
                end
                save([SavePath '/AutoSave_' 'AllElementsParameters.mat'],'AllElementsParameters');
                
                if LLcr_cats(j)/LLprev_cats(j)>0.95
                    StallInit2(j) = StallInit2(j)+1;
                end
                if (LLcr_cats(j)/LLprev_cats(j))>= StopCr && (OI  || StallInit2(j)<StallVal2(j))
                    stopped_learning_1(elem_j) = 1;
                    if ~isempty(Stored_QParam{model_i, elem_j})
                        AllElementsParameters{Index,2} = Stored_QParam{model_i, elem_j};
                    end
                    if IncludeStructuralAtt
                        if ~isempty(Stored_RegressionModel{model_i, elem_j})
                            AllElementsParameters{Index,4}  = Stored_RegressionModel{model_i, elem_j};
                        end
                    end
                    save([SavePath '/AutoSave_' 'learning_status_1.mat'],'stopped_learning_1');
                end
                if sum(stopped_learning_1) == total_num_elements
                    keep_learning = 0;
                end
                
                LLcr_cat = LLcr_cat + sum(LogLikVal);
                %% Testing with Testing set
                if IncludeStructuralAtt
                    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelTest,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,...
                        SpeedConstraints,1,ElementData.ModelTest.StrucAtt,KRparam(2:end),KRparam(1),...
                        KernelParameters{1},KernelParameters{2},[],[],[],AnnModel);
                    TestLogLik=sum(LogLikVal);
                else
                    [~,~,~,~,~,LogLikVal]=AnalysisObjective(ElementData.ModelTest,...
                        InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
                        Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                        OptimizationProceedure,OptLevel,...
                        SpeedConstraints,1);
                    TestLogLik=sum(LogLikVal);
                end
                if joint_inspector_train || stopped_learning_1(elem_j)
                    LLprev_cats(j) = LLcr_cats(j);
                    break
                end
            end
            j=j+1;
        end
        
    end
    LLcr = LLcr_cat;
    LLcr_cat = 0;
end
% check for place holders and reomve them
place_holder_id = find(strcmp(cellfun(@class,AllElementsParameters(:,2),'UniformOutput',false),'cell'));
if ~isempty(place_holder_id)
    for phi_ind = 1:length(place_holder_id)
        AllElementsParameters{place_holder_id(phi_ind),2}=[];
    end
end
