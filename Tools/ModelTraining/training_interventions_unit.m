d=uiprogressdlg(app.ModelTrainingUIFigure,'Title','Please Wait',...
            'Message',sprintf('Learning effect of interventions %d/%d',1,NumCats),'Cancelable','on'...
            ,'CancelText','Stop');
        
% time step
dt=1;

% transition for the kinimatic model
Ax1=[1 dt dt^2/2;
    0 1 dt;
    0 0 1];

A=blkdiag(Ax1,eye(3));

% Apply Constraints
ConstrainedKF = 1;

% Intervention Anlayese
InterventionCheck = 1;

% Observation matrix
F =[1 zeros(1,5)];
Material_index = ones(1,length(app.Tree.Children(2).Children));
for model_i=1:length(app.Tree.Children) % model type iterations: 2: SSM-KR or 1: SSM
    for elem_j=1:length(app.Tree.Children(model_i).Children)
        ElementName= app.Tree.Children(model_i).Children(elem_j).Text;
        d.Value = (elem_j/NumCats)/1.2;
        d.Message=sprintf('Strated learning effect of interventions on category %s',ElementName);
        disp("====================================================================================")
        disp(d.Message)
        if d.CancelRequested
            break
        end
        pause(0.001);
        % Read data
        app.Tree.SelectedNodes=app.Tree.Children(model_i).Children(elem_j);
        ElementData=load([FullPath 'TrainingData_Intervention_' erase(ElementName,"/") '.mat']);
        ElementData=struct2cell(ElementData);
        ElementData=ElementData{1};
        Index=find(strcmp(app.AllElementsParameters(:,1),ElementName));
        Ncurve = AllElementsParameters{Index,10};
        for j=1:size(ElementData.YS,3)
            if j==1
                [~,ElementData.init_x]=SpaceTransformationVec(Ncurve,ElementData.init_x,100,25);
                [~,ElementData.ModelValid.init_x]=SpaceTransformationVec(Ncurve,ElementData.ModelValid.init_x,100,25);
                [~,ElementData.ModelTest.init_x]=SpaceTransformationVec(Ncurve,ElementData.ModelTest.init_x,100,25);
            end
            [~,ElementData.YS(1,:,j)]=SpaceTransformationVec(Ncurve,ElementData.YS(1,:,j),100,25);
            [~,ElementData.ModelValid.YS(1,:,j)]=SpaceTransformationVec(Ncurve,ElementData.ModelValid.YS(1,:,j),100,25);
            [~,ElementData.ModelTest.YS(1,:,j)]=SpaceTransformationVec(Ncurve,ElementData.ModelTest.YS(1,:,j),100,25);
        end
        % remove unnecessary attributes
        if sum(strcmp(ElementName,{'Bande médiane', 'Tirants', 'Élément en élastomère','Toiture'})) == 1
            ElementData.StrucAtt(:,:,1) = [];
            ElementData.ModelValid.StrucAtt(:,:,1) = [];
            ElementData.ModelTest.StrucAtt(:,:,1) = [];
            Material_index(elem_j) = 0;
        end
        Yearly=ElementData.yearlyS;
        InspectorIDLabel=ElementData.InspectorLabelS;
        [~,MAxCondition]=SpaceTransformation(Ncurve,100,100,25);
        
        %% KR
        if model_i==2
            AttStruc=reshape((ElementData.StrucAtt),size(ElementData.StrucAtt,2),size(ElementData.StrucAtt,3));
            AvgObs=reshape((ElementData.init_x),size(ElementData.init_x,2),1);
            AllAtt=[AttStruc AvgObs];
            RegressionModel = AllElementsParameters{Index,4};
            %RegressionModel.categ_st_att_ind = Material_index(elem_j);
        else
            RegressionModel = [];
            AllAtt = [];
        end
        
        param = AllElementsParameters{Index,2};
        
        Q_ki = param(1).^2*[(dt^4)/4,(dt^3)/2,(dt^2)/2;
                        (dt^3)/2,(dt^2)/1,(dt^1)/1;
                        (dt^2)/2,(dt^1)/1,1];
        Q=blkdiag(Q_ki,zeros(3));
        
        % inspection data
        y = (ElementData.YS);
        
        % Noise
        Re = ElementData.ReS;
        Be = ElementData.InpecBiaseS;

        % Intervention Vector
        InterventionVector = (ElementData.InterventionVector);
        
        init_x=zeros(6,1);

        IntervType = (ElementData.InterventionType);
        IntervType = fix(IntervType./10.^fix(log10(IntervType)));
        
        InspectorsData = AllElementsParameters(Index,3);
        Ncurve = AllElementsParameters{Index,10};
        OptBounds_int = {OptBoundsInt(1,:) OptBoundsInt(2,:) OptBoundsInt(3,:) ...
            OptBoundsInt(4,:) OptBoundsInt(5,:) OptBoundsInt(6,:)};
        ModelParamLocal=app.AllElementsParameters{Index,2}{8}(:,1)';
        for IndexTypeVal=1:3
            Cindex = find(IntervType==IndexTypeVal);
            if ~isempty(Cindex)
                InspectorsID=(ElementData.InspectorLabelS);
                [ModelParamLocal,~,~,fx_NR]=Newton_Raphson(@(ModelParamLocal)...
                    opt_intervention(Cindex,InspectorsID,InspectorsData,Re,Be,y,...
                    param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
                    InterventionVector,model_i,RegressionModel,AllAtt,IndexTypeVal)...
                    ,ModelParamLocal,'log_transform','no','output','original','laplace','no',...
                    'convergence_tol',1E-1,'bounds',OptBounds_int);
                [~,InterventionMu_Network,InterventionVar_Network]= ...
                    opt_intervention(Cindex,InspectorsID,InspectorsData,Re,Be,y,...
                    param,ModelParamLocal,A,F,Q,Ncurve,ConstrainedKF,InterventionCheck,...
                    InterventionVector,model_i,RegressionModel,AllAtt,IndexTypeVal);
                int_Ex{1,IndexTypeVal} = InterventionMu_Network;
                int_Ex{2,IndexTypeVal} = InterventionVar_Network;
                int_param{1,IndexTypeVal} = ModelParamLocal;
            else
                int_Ex{1,IndexTypeVal} = nan(3,1);
                int_Ex{2,IndexTypeVal} = nan(3);
                int_param{1,IndexTypeVal} = nan(1,6);
            end
            AllElementsParameters{Index,6}{1,IndexTypeVal}=int_Ex{1,IndexTypeVal};
            AllElementsParameters{Index,7}{1,IndexTypeVal}=int_Ex{2,IndexTypeVal};
            AllElementsParameters{Index,8}{1,IndexTypeVal}=int_param{1,IndexTypeVal};
        end        
    end
end