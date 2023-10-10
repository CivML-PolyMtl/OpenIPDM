d=uiprogressdlg(app.ModelTrainingUIFigure,'Title','Please Wait',...
            'Message',sprintf('Learning parameters of structural category %d/%d',1,NumCats),'Cancelable','on'...
            ,'CancelText','Stop');

InspStrucIndex=[];

%% Kalman Filter
dt=1;
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
% Q=@(param) param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;
%                         (dt^4)/8 (dt^3)/3 (dt^2)/2;
%                         (dt^3)/6 (dt^2)/2 dt];
                    
Q=@(param) param(1).^2*[(dt^4)/4,(dt^3)/2,(dt^2)/2;
                        (dt^3)/2,(dt^2)/1,(dt^1)/1;
                        (dt^2)/2,(dt^1)/1,1];

F=[1,0,0]; 
x0=zeros(3,1);
s2_X0=zeros(3);
SpeedConstraints=[1 -50 0];
InitN=OptBoundsData(1,2);
FinN=OptBoundsData(1,3);
LogLikCR=-10^10;
cla(app.UIAxes);
if ResumeTraining
    InitN = AllElementsParameters{Index,10};
end
LL_validation = 0;
LL_validation_prev = -10^8;

AllElementsParameters = app.AllElementsParameters;
Material_index = ones(1,length(app.Tree.Children(2).Children));
for Ncurve=InitN:FinN
    for model_i=1:length(app.Tree.Children) % model type iterations: 2: SSM-KR or 1: SSM
        for elem_j=1:length(app.Tree.Children(model_i).Children)
            ElementName= app.Tree.Children(model_i).Children(elem_j).Text;
            d.Value = (elem_j/NumCats)/1.2;
            d.Message=sprintf('Started training structural category %s',ElementName);
            disp(d.Message)
            if d.CancelRequested
                break
            end
            pause(0.001);
            % Read data
            
            app.Tree.SelectedNodes=app.Tree.Children(model_i).Children(elem_j);
            ElementData=load([FullPath 'TrainingData_' erase(ElementName,"/") '.mat']);
            ElementData=struct2cell(ElementData);
            ElementData=ElementData{1};
            Index=find(strcmp(app.AllElementsParameters(:,1),ElementName));
            OptBoundsData=app.AllElementsParameters{Index,2}{1};
            % Inspector uncertainty bounds
            bounds={OptBoundsData(3,2:3) [-OptBoundsData(3,3) OptBoundsData(3,3)]};

            % Model parameter bounds
            boundsQ={OptBoundsData(2,2:3) OptBoundsData(3,2:3) OptBoundsData(4,2:3)...
                OptBoundsData(5,2:3) OptBoundsData(6,2:3) OptBoundsData(7,2:3)};
            InspectorsData=load([FullPath 'Inspectors_' erase(ElementName,"/") '.mat']);
            InspectorsData=struct2cell(InspectorsData);
            InspectorsData{1}(:,2)= 0 ;
            InspectorsData{1}(:,3)= OptBoundsData(3,1);
            AllElementsParameters{Index,3} = InspectorsData{1};
            config_opt.elements = app.AllElementsParameters;
            config_opt.indexes = [model_i;elem_j];
            save([SavePath '/AutoSave_Config_' erase(ElementName,"/") '.mat'],'config_opt');
            % Initial optimization step
            OptLevel=1;
            OptimizationProceedure=1;
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
            elseif strcmp(app.Switch.Value, 'Synthetic')
                Material_index(elem_j) = 0;
            end
            save([FullPath 'TrainingData_' erase(ElementName,"/") '_' sprintf('%d',Ncurve) '.mat'],'ElementData', '-v7.3');
            if ResumeTraining
                PARAM = AllElementsParameters{Index,2};
                InspectorsData{1} = AllElementsParameters{Index,3};
                Stored_InspectorData = InspectorsData;
                Stored_QParam = PARAM;
            else
                % initial parameter values
                PARAM=[OptBoundsData(2,1) OptBoundsData(3,1) OptBoundsData(4,1) ...
                    OptBoundsData(5,1) OptBoundsData(6,1) OptBoundsData(7,1)];
                [PARAM,~,~,fx_NR]=Newton_Raphson(@(PARAM) AnalysisObjective(ElementData,...
                    InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
                    x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
                    OptimizationProceedure,OptLevel,SpeedConstraints,1)...
                    ,PARAM,'log_transform','no','output','original','laplace','no',...
                    'convergence_tol',1E-1,'bounds',boundsQ);
%                 PARAM = [0.0055675      2.3604      2.5032    0.049997    0.025623    0.060676];
%                 fx_NR = -10e6;
                AllElementsParameters{Index,3}(:,3)=PARAM(2);
                AllElementsParameters{Index,3}(:,2)=0;
                LogLik=fx_NR(end);%-10000;%
                Stored_QParam = PARAM;
                Stored_InspectorData = [];
            end
            % save parameters
            AllElementsParameters{Index,2}=PARAM;
            AllElementsParameters{Index,10}=Ncurve;
            % reset the Ann for the next transformation parameter
            if Ncurve > InitN
                AnnModel=[];
                NumAttributes=0;
                if model_i==1
                    IncludeStructuralAtt=0;
                else
                    IncludeStructuralAtt=1;
                    if strcmp(app.Switch.Value, 'Synthetic')
                        NumAttributes = 1;
                    elseif sum(strcmp(ElementName,{'Chasse-roue / trottoir','Tirants','Revêtement de mur','Voûte / Dalle','Toiture'})) == 1
                        NumAttributes = 10;
                    elseif sum(strcmp(ElementName, {'Dessous de la dalle/voûte','Mur de tête','Mur en aile','Murs/naiss.voûte/coins infér.'}))== 1
                        NumAttributes = 9;
                    elseif strcmp(ElementName,'Mur')
                        NumAttributes = 6;
                    else 
                        NumAttributes = 11;
                    end
                    app.AllElementsParameters{Index,2}{5} = NumAttributes;
                    % TAGI parameters
                    StructuralAttributes=0;
                    if isempty(AnnModel)
                        AnnModel.batch_size = 2;%TAGI_params{2};
                        AnnModel.max_epoches = 30;%TAGI_params{1};
                        AnnModel.cont_training = 0;
                        AnnModel.evaluate_mode = 1;
                        AnnModel.theta_nn = [];
                    end
                end
            end
            AllElementsParameters{Index,4}=AnnModel;
        end
    end
    save([SavePath '/AutoSave_' 'AllElementsParameters.mat'],'AllElementsParameters');
    Stage2_unit();
    LL_validation = sum(LLprev_cats);
    if LL_validation > LL_validation_prev
        store_AllElementsParameters = AllElementsParameters;
    else
        AllElementsParameters = store_AllElementsParameters;
    end

end



