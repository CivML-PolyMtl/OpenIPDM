% Intervention Model parameter bounds
bounds={OptBoundsInt(1,2:3) OptBoundsInt(2,2:3) OptBoundsInt(3,2:3)...
    OptBoundsInt(4,2:3) OptBoundsInt(5,2:3) OptBoundsInt(6,2:3)};


InspStrucIndex=[];

%% Kalman Filter
dt=1;
A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
Q_ki=@(PARAMstore) PARAMstore(1).^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;
                        (dt^4)/8 (dt^3)/3 (dt^2)/2;
                        (dt^3)/6 (dt^2)/2 dt];
Q=blkdiag(Q_ki,zeros(3));
InterventionVar_Network=diag([ModelParamLocal(4)^2 ModelParamLocal(5)^2 ModelParamLocal(6)^2]);
Q_r=diag([ModelParamLocal(1)^2 ModelParamLocal(2)^2 ModelParamLocal(3)^2]);
Q_r=blkdiag(Q_r,Q_r);
F=[1,zeros(1,5)]; 
x0=zeros(6,1);
s2_X0=zeros(6);
SpeedConstraints=[1 -50 0];
LogLikCR=-10^10;

% Initial optimization step
OptimizationProceedure=1;
for j=1:size(ElementData.YS,3)
    if j==1
        [~,ElementData.init_x]=SpaceTransformation(Ncurve,ElementData.init_x,100,25);
    end
    [~,ElementData.YS(1,:,j)]=SpaceTransformation(Ncurve,ElementData.YS(1,:,j),100,25);
end
if IncludeStructuralAtt
    OptLevel=3;
    boundKr=[];
    InitialEx=[];
    InitialVar=[];
    Kr_Param=KernelParameters{3}(:,1)';
    KRparam=zeros(2,1);
    StructuralAttributes=0;
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
end
[Qparam,~,~,LogLik]=Newton_Raphson(@(PARAM) AnalysisObjective(ElementData,...
    InspectorsData{1}(:,1),InspStrucIndex,OptBoundsData,A,F,...
    Q, x0, s2_X0,PARAM,[0 0],InspectorsData{1},[],Ncurve,...
    OptimizationProceedure,OptLevel,SpeedConstraints,1,...
    StructuralAttributes,KRparam(2:end),KRparam(1),...
    KernelParameters{1},KernelParameters{2},InitialEx,InitialVar,[]),...
    PARAM,'log_transform','no','output','original','laplace',...
    'no','convergence_tol',1E-4,'bounds',boundsQ);



    
    
    
    


