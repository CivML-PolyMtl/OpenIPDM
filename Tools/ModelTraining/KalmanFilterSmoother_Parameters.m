init_x(1,:)=gather(y(1,:,2));
OriginalValues=[25:1:100];
[~,MAxCondition]=SpaceTransformation(Ncurve,OriginalValues,100,25);
NormObsValues=init_x(1,:);
NormObsValues(find(NormObsValues>MAxCondition(end)))=MAxCondition(end);
DifferenceObs=MAxCondition(end)-NormObsValues;
% init_x(1,:)=IniObsAvg(1,:);%gather(y(1,:,2));
if OptProcedure==1
    init_x(2,:)=0;%param(6)*(DifferenceObs);%+param(8);
    init_x(3,:)=0;%param(7)*init_x(2,:);
else
    init_x(2,:)=-0.01;
    init_x(3,:)=-0.001;
end
if OptProcedure==1
    if GPUCompute
        init_V=zeros(3,3,length(init_x(1,:)),'gpuArray');
    else
        init_V=zeros(3,3,1);
    end
    if OptLevel~=1
        param(3)=param(2);
    end
    init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
    init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);%Re(1,:,2).*param(6).^2+param(4)^2;%
    init_V(3,3,:)=param(5).^2;%(param(5).^2).*init_V(2,2,:)+init_V(2,2,:).*param(7)^2+
    
end
[x, Var, ~, ~,~,~] = kalman_filter(y, A, C, Q, R, Re...
    , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,ObsYears...
    ,Ncurve,OptBoundsData,SpeedConstraints(1,1),SpeedConstraints,GPUCompute,Nsigma);

TotalTimeSteps=size(y,3);
TrainSmootherRun();
if GPUCompute
    InitialCond=gather(ExSmooth(1,:,2));
    InitialCond(find(InitialCond>MAxCondition(end)))=MAxCondition(end);
    DifferenceObs=MAxCondition(end)-InitialCond;
else
    InitialCond=ExSmooth(1,2);
    InitialCond(find(InitialCond>MAxCondition(end)))=MAxCondition(end);
    DifferenceObs=MAxCondition(end)-InitialCond;
end
if OptProcedure==1
    init_x(2,:)=0;%param(6)*(DifferenceObs);%+param(8);
    init_x(3,:)=0;%param(7)*init_x(2,:);
else
    init_x(2,:)=-0.01;
    init_x(3,:)=-0.001;
end
if OptProcedure==1
    if GPUCompute
        init_V=zeros(3,3,length(init_x(1,:)),'gpuArray');
    else
        init_V=zeros(3,3,1);
    end
    init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
    init_V(2,2,:)=(param(4).^2).*(DifferenceObs)+(param(6).^2);%reshape(gather(VarSmooth(1,1,:,2)),length(VarSmooth(1,1,:,2)),1).*param(6).^2;%
    init_V(3,3,:)=(param(5).^2);%.*init_V(2,2,:)+;init_V(2,2,:).*param(7)^2;%
end
[~, ~, ~, loglik,~,~] = kalman_filter(y, A, C, Q, R, Re...
    , init_x, init_V,InpecBiase,CurrentInspectorObs,RU,InspBU,...
    ObsYears,Ncurve,OptBoundsData,SpeedConstraints(1,1),...
    SpeedConstraints,GPUCompute,Nsigma);