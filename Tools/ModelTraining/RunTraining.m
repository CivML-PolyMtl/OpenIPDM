
% initial parameter values
PARAM=[OptBoundsData(2,1) OptBoundsData(3,1) OptBoundsData(4,1) ...
    OptBoundsData(5,1) OptBoundsData(6,1) OptBoundsData(7,1)];

% Inspector uncertainty bounds
bounds={OptBoundsData(3,2:3) [-OptBoundsData(3,3) OptBoundsData(3,3)]};

% Model parameter bounds
boundsQ={OptBoundsData(2,2:3) OptBoundsData(3,2:3) OptBoundsData(4,2:3)...
    OptBoundsData(5,2:3) OptBoundsData(6,2:3) OptBoundsData(7,2:3)};


InspStrucIndex=[];

%% Kalman Filter
dt=1;
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q=@(param) param(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;
                        (dt^4)/8 (dt^3)/3 (dt^2)/2;
                        (dt^3)/6 (dt^2)/2 dt];
F=[1,0,0]; 
x0=zeros(3,1);
s2_X0=zeros(3);
SpeedConstraints=[1 -50 0];
InspectorsData{1}(:,2)=OptBoundsData(3,1);
InitN=OptBoundsData(1,2);
FinN=OptBoundsData(1,3);
LogLikCR=-10^10;
cla(app.UIAxes);
if ResumeTraining
    InitN = AllElementsParameters{Index,10};
end
for Ncurve=InitN:FinN
    % Read data
    ElementData=load([FullPath 'TrainingData_' ElementName '.mat']);
    ElementData=struct2cell(ElementData);
    ElementData=ElementData{1};
    % Initial optimization step
    OptLevel=1;
    OptimizationProceedure=1;
    for j=1:size(ElementData.YS,3)
        if j==1
            [~,ElementData.init_x]=SpaceTransformation(Ncurve,ElementData.init_x,100,25);
        end
        [~,ElementData.YS(1,:,j)]=SpaceTransformation(Ncurve,ElementData.YS(1,:,j),100,25);
    end
    if ResumeTraining
        PARAM = AllElementsParameters{Index,2};
        InspectorsData{1} = AllElementsParameters{Index,3};
    else
        [PARAM,~,~,fx_NR]=Newton_Raphson(@(PARAM) AnalysisObjective(ElementData,...
            InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F, Q,...
            x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
            OptimizationProceedure,OptLevel,SpeedConstraints,1)...
            ,PARAM,'log_transform','no','output','original','laplace','no',...
            'convergence_tol',1E-3,'bounds',boundsQ);
        InspectorsData{1}(:,3)=PARAM(2);
        InspectorsData{1}(:,2)=0;
        LogLik=fx_NR(end);% LogLik=0;%
    end
    
    
    
    if IncludeStructuralAtt
        if ResumeTraining 
            InitialEx = AllElementsParameters{Index,4}.InirilizedEx;
            InitialVar = AllElementsParameters{Index,4}.InirilizedVar;
            Kr_Param = [AllElementsParameters{Index,4}.Sigma_W0... 
                        AllElementsParameters{Index,4}.Kernel_l];
            StructuralAttributes=ElementData.StrucAtt;
            KRparam = Kr_Param';
        else
            InitialEx=[];
            InitialVar=[];
            Kr_Param=KernelParameters{3}(:,1)';
            StructuralAttributes=0;
            KRparam=zeros(2,1);
        end
        OptLevel=3;
        boundKr=[];
        
        
        for vi=1:length(KernelParameters{3}(:,1))
            boundKr=[boundKr {KernelParameters{3}(vi,2:3)}];
        end
    end
    
    % Factoring Structural Attributes
    % IncludeStructuralAtt
    if ~IncludeStructuralAtt
        KRparam=zeros(2,1);
        KernelParameters{1}=0;
        KernelParameters{2}=0;
        InitialEx=[];
        InitialVar=[];
        RegressionModel = [];
    end
    
    % Full optimization
    Stage2();
    if TestLogLik>LogLikCR
        LogLikCR=TestLogLik;
        Nstore=Ncurve;
        if IncludeStructuralAtt
            RegressionModelstore=RegressionModel;
        else
            RegressionModelstore=[];
        end
        PARAMstore=PARAM;
    end
end

