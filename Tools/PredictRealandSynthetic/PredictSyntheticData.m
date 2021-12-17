function PredictSyntheticData(SynDatabaseShort,SynDatabaseState,M,...
    SynInspTrue,SynInsp,ObsWindow,NumOfElements,RegressionModel,TableOfParameters,...
    GraphMatrix,NumOfFigures)
if isempty(GraphMatrix)
    GraphMatrix=zeros(5);
end
LogLik=0;
Totallength=1+ObsWindow+M;
BoundViolationAll=NaN(NumOfElements,Totallength);
dt=1;
dfct_n=@(xts,nts)(nts*exp(-xts.^nts))/gamma(1/nts);
if NumOfFigures>0
    RndID=randi([1,NumOfElements],NumOfFigures,1);
end
if iscell(SynInsp)
    SynInsp=cell2mat(SynInsp);
    SynInsp=fliplr(SynInsp);
    %     SynInsp(:,2)=[];
end
for ll=1:NumOfElements
    % Element ID
    ElementInd=ll;
    param=TableOfParameters{2,1};
    paramTrue=TableOfParameters{1,1};
    
    A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
    Q= param.^2*[(dt^4)/4 (dt^3)/2 (dt^2)/2;(dt^3)/2 (dt^2)/1 (dt^1)/1;(dt^2)/2 (dt^1)/1 1];
    QTrue= paramTrue.^2*[(dt^4)/4 (dt^3)/2 (dt^2)/2;(dt^3)/2 (dt^2)/1 (dt^1)/1;(dt^2)/2 (dt^1)/1 1];
    C=[1,0,0];
    % initial knowledge at time t=0
    init_x=[100; TableOfParameters{2,5}; TableOfParameters{2,6}];
    init_V=[TableOfParameters{2,2}^2,0,0;0,TableOfParameters{2,3},0;0,0,TableOfParameters{2,4}];
    
    y=SynDatabaseShort{ElementInd}(:,1);
    StructureInd=SynDatabaseShort{ElementInd}(1,4);
    yearly=SynDatabaseShort{ElementInd}(1,3):SynDatabaseShort{ElementInd}(end,3)+M;
    [~,~,iinsp]=intersect(SynDatabaseShort{ElementInd}(:,3),yearly);
    ObsYears=zeros(length(yearly),1);
    ObsYears(iinsp)=1;
    Diffy=nan(length(yearly),1);
    Diffy(iinsp)=y;
    yinsp=y;
    y=Diffy';
    ReTrue=[];
    % 0 is any inspector
    % NaN is Unkown inspector (prediction)
    Insp=SynDatabaseShort{ElementInd}(:,2)';
    for k=1:length(Insp)
        [~,ia,~] = intersect(SynInsp(:,end),Insp(k));
%         if SynInsp(ia,3) == 99173
%             ll;
%         end
        Re(k)=SynInsp(ia,1).^2;
        InpecBiase(k)=SynInsp(ia,2);
        ReTrue(k)=SynInspTrue(ia,1).^2;
        InpecBiaseTrue(k)=SynInspTrue(ia,2);
    end
    ReZero=zeros(1,length(y));
    ReZeroTrue=zeros(1,length(y));
    InpecBiaseZero=zeros(1,length(y));
    InpecBiaseTrueZero=zeros(1,length(y));
    InspectorLabel=nan(1,length(y));
    ReZero(iinsp)=Re;
    InpecBiaseZero(iinsp)=InpecBiase;
    InspectorLabel(iinsp)=Insp;
    Re=ReZero;
    ReZeroTrue(iinsp)=ReTrue;
    ReTrue=ReZeroTrue;
    InpecBiase=InpecBiaseZero;
    InpecBiaseTrueZero(iinsp)=InpecBiaseTrue;
    InpecBiaseTrue = InpecBiaseTrueZero;
    OptmInsp=zeros(1,length(y));
    RU=1;InspBU=0;R=0;
    if sum(find(ll==RndID))
        FigureID=1;
    else
        FigureID=0;
    end
    if ~isempty(RegressionModel)
        ObsIndexes=find(~isnan(y));                                                 % Identify the indecies of observations (to take the average)
        if length(ObsIndexes)==2
            ObsIndexes=ObsIndexes(1:2);
        elseif length(ObsIndexes)>2
            ObsIndexes=ObsIndexes(1:3);
        else
            ObsIndexes=ObsIndexes(1);
        end
        AllAtt=[SynDatabaseShort{ElementInd}(1,end) mean(y(ObsIndexes))];
    else
        AllAtt=[];
    end
    [loglik,MAcc,AccVal,AccST,PWithinCI,PWithinCITrue,Ebar,Ebarbar,...
        BoundViolation,XtrueValue,XValue,XbtrueTr,XbbtrueTr,XestimateTr,~,~,...
        Erbias, Erbarbias, AccValBias]=KFsKF([],y, A, C, Q, R, Re,init_x, init_V,InpecBiase,InpecBiaseTrue,OptmInsp,RU,...
        InspBU,ObsYears,yearly,InspectorLabel,StructureInd,ElementInd,...
        TableOfParameters{2,8},ReTrue,QTrue,SynDatabaseState,RegressionModel,AllAtt,...
        TableOfParameters,GraphMatrix,[],FigureID);
    LogLik=LogLik+loglik;
    LogLikVec(ll)=loglik;
    XtrueValueAll{1,ll}=XtrueValue;
    XValueAll{1,ll}=XValue(1,:);
    EbarAll{1,ll}=Ebar;
    EbarbarAll{1,ll}=Ebarbar;
    AccAll{1,ll}=AccVal;
    ErbiasAll{1,ll}=Erbias;
    ErbarbiasAll{1,ll}=Erbarbias;
    AccValBiasAll{1,ll}=AccValBias;
    AccAllSTD{1,ll}=AccST;
    PAllWithinCI{1,ll}=PWithinCI;
    PAllWithinCITrue{1,ll}=PWithinCITrue;
    SIndex=length(AccAll{1,ll})-M+1;
    S2Index=Totallength-(length(AccAll{1,ll}))+1;
    BoundViolationAll(ll,S2Index:length(BoundViolation)+S2Index-1)=BoundViolation;
    xTr(ll)=interp1([-12.5,137.5],[-2,2],XestimateTr(1,SIndex-1));
    XbValueAll(1,ll)=XestimateTr(2,SIndex-1).*dfct_n(xTr(ll),2^TableOfParameters{2,8});
    XbbValueAll(1,ll)=XestimateTr(3,SIndex-1);
    XbtrueTrAll(1,ll)=XbtrueTr(SIndex-1);
    XbbtrueTrAll(1,ll)=XbbtrueTr(SIndex-1);
    % initial speed and acc
    xTr1(ll)=interp1([-12.5,137.5],[-2,2],XestimateTr(1,1));
    XbValueAllinit(1,ll)=XestimateTr(2,1).*dfct_n(xTr1(ll),2^TableOfParameters{2,8});
    XbbValueAllinit(1,ll)=XestimateTr(3,1);
    XbtrueTrAllinit(1,ll)=XbtrueTr(1);
    XbbtrueTrAllinit(1,ll)=XbbtrueTr(1);
    
    AccAll{1,ll}=AccAll{1,ll}(SIndex:end);
    XtrueValueAll{1,ll}=XtrueValueAll{1,ll}(SIndex:end);
    XValueAll{1,ll}=XValueAll{1,ll}(SIndex:end);
    AccAllSTD{1,ll}=AccAllSTD{1,ll}(SIndex:end);
    PAllWithinCI{1,ll}=PAllWithinCI{1,ll}(SIndex:end);
    PAllWithinCITrue{1,ll}=PAllWithinCITrue{1,ll}(SIndex:end);
    EbarAll{1,ll}=EbarAll{1,ll}(SIndex:end);
    EbarbarAll{1,ll}=EbarbarAll{1,ll}(SIndex:end);
    ErbiasAll{1,ll}=ErbiasAll{1,ll}(SIndex:end);
    ErbarbiasAll{1,ll}=ErbarbiasAll{1,ll}(SIndex:end);
    AccValBiasAll{1,ll}=AccValBiasAll{1,ll}(SIndex:end);
    BoundViolationAllTime=nanmean(BoundViolationAll);
    for vi=1:length(AccAll{1,ll})
        AccAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), AccAll(1,:)));
        AccAllTime_std(vi,1)=std(cellfun(@(c) c(1,vi), AccAll(1,:)));
        AccAllSTDTime(vi,1)=mean(cellfun(@(c) c(1,vi), AccAllSTD(1,:)));
        PAllWithinCITime(vi,1)=mean(cellfun(@(c) c(1,vi), PAllWithinCI(1,:)));
        PAllWithinCITrueTime(vi,1)=mean(cellfun(@(c) c(1,vi), PAllWithinCITrue(1,:)));
        EbarAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), EbarAll(1,:)));
        EbarAllTime_std(vi,1)=mean(cellfun(@(c) c(1,vi), EbarAll(1,:)));
        EbarbarAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), EbarbarAll(1,:)));
        EbarbarAllTime_std(vi,1)=mean(cellfun(@(c) c(1,vi), EbarbarAll(1,:)));
        ErbiasAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), ErbiasAll(1,:)));
        ErbarbiasAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), ErbarbiasAll(1,:)));
        ErbiasAllTime_Std(vi,1)=std(cellfun(@(c) c(1,vi), ErbiasAll(1,:)));
        ErbarbiasAllTime_Std(vi,1)=std(cellfun(@(c) c(1,vi), ErbarbiasAll(1,:)));
        AccValBiasAllTime(vi,1)=mean(cellfun(@(c) c(1,vi), AccValBiasAll(1,:)));
        AccValBiasAllTime_Std(vi,1)=std(cellfun(@(c) c(1,vi), AccValBiasAll(1,:)));
    end
    clear y Re ReZero StructureInd yearly Insp RU InspBU R InpecBiase InpecBiaseTrue
end
if GraphMatrix(5,2)
    figure(1)
    hold off
    plot(1:length(BoundViolationAllTime),BoundViolationAllTime);
    xlabel('Time (Years)')
    ylabel('CI with Increasing Pattern')
    ylim([0,1])
    grid on
end

if GraphMatrix(2,1)
    figure(2)
    Tval=1:length(EbarbarAllTime);
    subplot(1,3,1)
    hold off
    plot(1:length(EbarbarAllTime),EbarbarAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[EbarbarAllTime'+EbarbarAllTime_std',fliplr(EbarbarAllTime'-EbarbarAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[EbarbarAllTime'+2*EbarbarAllTime_std',fliplr(EbarbarAllTime'-2*EbarbarAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $|\ddot{x}_{t,p}^j -\ddot{\mu}_{t|\mathtt{T},p}^j|$','interpreter','latex')
    ylim([0,0.1])
    grid on
    
    subplot(1,3,2)
    hold off
    plot(1:length(EbarAllTime),EbarAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[EbarAllTime'+EbarAllTime_std',fliplr(EbarAllTime'-EbarAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[EbarAllTime'+2*EbarAllTime_std',fliplr(EbarAllTime'-2*EbarAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $|\dot{x}_{t,p}^j -\dot{\mu}_{t|\mathtt{T},p}^j|$','interpreter','latex')
    ylim([0,3])
    grid on
    
    subplot(1,3,3)
    hold off
    plot(1:length(AccAllTime),AccAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[AccAllTime'+AccAllTime_std',fliplr(AccAllTime'-AccAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[AccAllTime'+2*AccAllTime_std',fliplr(AccAllTime'-2*AccAllTime_std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $|{x}_{t,p}^j -{\mu}_{t|\mathtt{T},p}^j|$','interpreter','latex')
    grid on
    
    figure(3)
    Tval=1:length(AccValBiasAllTime);
    subplot(1,3,1)
    hold off
    plot(1:length(ErbarbiasAllTime),ErbarbiasAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[ErbarbiasAllTime'+ErbarbiasAllTime_Std',fliplr(ErbarbiasAllTime'-ErbarbiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[ErbarbiasAllTime'+2*ErbarbiasAllTime_Std',fliplr(ErbarbiasAllTime'-2*ErbarbiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $(\ddot{x}_{t,p}^j -\ddot{\mu}_{t|\mathtt{T},p}^j)$','interpreter','latex')
    ylim([-0.2,0.2])
    grid on
    
    subplot(1,3,2)
    hold off
    plot(1:length(ErbiasAllTime),ErbiasAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[ErbiasAllTime'+ErbiasAllTime_Std',fliplr(ErbiasAllTime'-ErbiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[ErbiasAllTime'+2*ErbiasAllTime_Std',fliplr(ErbiasAllTime'-2*ErbiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $(\dot{x}_{t,p}^j -\dot{\mu}_{t|\mathtt{T},p}^j)$','interpreter','latex')
    ylim([-2,2])
    grid on
    
    subplot(1,3,3)
    hold off
    plot(1:length(AccValBiasAllTime),AccValBiasAllTime);
    hold on
    patch([Tval,fliplr(Tval)],[AccValBiasAllTime'+AccValBiasAllTime_Std',fliplr(AccValBiasAllTime'-AccValBiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    patch([Tval,fliplr(Tval)],[AccValBiasAllTime'+2*AccValBiasAllTime_Std',fliplr(AccValBiasAllTime'-2*AccValBiasAllTime_Std')],'k','FaceAlpha',0.2,'EdgeColor','none')
    xlabel('Time (Years)')
    ylabel('Average $({x}_{t,p}^j -{\mu}_{t|\mathtt{T},p}^j)$','interpreter','latex')
    ylim([-10,10])
    grid on
end

if GraphMatrix(4,2)
    figure(4)
    hold off
    plot(1:length(AccAllSTDTime),AccAllSTDTime);
    xlabel('Time (Years)')
    ylabel('Average $|x_t - \mu_{t|T}|/\sigma$','interpreter','latex')
    ylim([0,6])
    grid on
end

if GraphMatrix(3,2)
    figure(5)
    subplot(2,1,1)
    hold off
    plot(1:length(PAllWithinCITime),PAllWithinCITime);
    xlabel('Forecast Time (Years)')
    ylabel('Average $\Pr(|x^j_{t,p}-\mu^j_{t|\mathtt{T},p}|\le2\sigma^j_{t|\mathtt{T},p}$)','interpreter','latex')
    ylim([0,1])
    grid on
    
    subplot(2,1,2)
    hold off
    plot(1:length(PAllWithinCITrueTime),PAllWithinCITrueTime);
    xlabel('Forecast Time (Years)')
    ylabel('Average $\Pr(|x^j_{t,p}-\mu^j_{t|\mathtt{T},p}|\le2\sigma^j_{t|\mathtt{T},p}$)(true)','interpreter','latex')
    ylim([0,1])
    grid on
end

if GraphMatrix(3,1)
    figure(6)
    subplot(1,3,1)
    hold off
    plot(cellfun(@(c) c(1,1), XtrueValueAll(1,:)),cellfun(@(c) c(1,1),...
        XValueAll(1,:)),'.');
    xlabel('True (1 Year)')
    ylabel('Estimated (1 Year)')
    xlim([25 100])
    ylim([25 100])
    axis square
    grid on
    
    subplot(1,3,2)
    hold off
    plot(cellfun(@(c) c(1,5), XtrueValueAll(1,:)),cellfun(@(c) c(1,5),...
        XValueAll(1,:)),'.');
    xlabel('True (5 Year)')
    ylabel('Estimated (5 Year)')
    xlim([25 100])
    ylim([25 100])
    axis square
    grid on
    
    subplot(1,3,3)
    hold off
    plot(cellfun(@(c) c(1,10), XtrueValueAll(1,:)),cellfun(@(c) c(1,10),...
        XValueAll(1,:)),'.');
    xlabel('True (10 Year)')
    ylabel('Estimated (10 Year)')
    xlim([25 100])
    ylim([25 100])
    axis square
    grid on
end

if GraphMatrix(4,1)
    figure(7)
    subplot(1,2,1)
    hold off
    plot(XbtrueTrAll(1,:),XbValueAll(1,:),'.');
    xlabel('$\dot{x}_{\mathtt{T}}$','interpreter','latex')
    ylabel('$\dot{\mu}_{\mathtt{T}|\mathtt{T}}$','interpreter','latex')
    xlim([-1 0])
    ylim([-1 0])
    hold on
    plot([-3,0],[-3,0]);
    axis square
    grid on
    
    figure(7)
    subplot(1,2,2)
    hold off
    plot(XbbtrueTrAll(1,:),XbbValueAll(1,:),'.');
    xlabel('$\ddot{x}_{\mathtt{T}}$','interpreter','latex')
    ylabel('$\ddot{\mu}_{\mathtt{T}|\mathtt{T}}$','interpreter','latex')
    xlim([-0.1 0])
    ylim([-0.1 0])
    hold on
    plot([-3,0],[-3,0]);
    axis square
    grid on
end

if GraphMatrix(2,2)
    figure(8)
    hold off
    histogram((XbtrueTrAll(1,:)-XbValueAll(1,:)));
    xlabel('Histogram $(\dot{x}_{\mathtt{T}}-\dot{\mu}_{\mathtt{T}|\mathtt{T}})$','interpreter','latex')
    grid on
end

if GraphMatrix(5,1)
    figure(9)
    subplot(1,2,1)
    hold off
    plot(XbtrueTrAllinit(1,:),XbValueAllinit(1,:),'.');
    xlabel('$\dot{x}_{0}$','interpreter','latex')
    ylabel('$\dot{\mu}_{0}$','interpreter','latex')
    xlim([-1 0])
    ylim([-1 0])
    hold on
    plot([-3,0],[-3,0]);
    axis square
    grid on
    
    figure(9)
    subplot(1,2,2)
    hold off
    plot(XbbtrueTrAllinit(1,:),XbbValueAllinit(1,:),'.');
    xlabel('$\ddot{x}_{0}$','interpreter','latex')
    ylabel('$\ddot{\mu}_{0}$','interpreter','latex')
    xlim([-0.1 0])
    ylim([-0.1 0])
    hold on
    plot([-3,0],[-3,0]);
    axis square
    grid on
end