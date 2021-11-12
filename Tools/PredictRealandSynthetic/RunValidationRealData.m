
if isempty(AllParameters)
    msgbox('Please Specify Model Parameters before Running the Analysis')
else
    ActiveParametersSet=AllParameters.ActiveParametersSet;
    TableEstimatedParameters=AllParameters.TableEstimatedParameters;
    RealInspectorsEstimatedSigma=AllParameters.RealInspectorsEstimatedSigma;
    RealInspectorsID=AllParameters.RealInspectorsID;
    if ~isempty(AllParameters.RegressionModel)
        RegressionModel=AllParameters.RegressionModel;
    else
        RegressionModel=[];
    end
    
    MdataEngyUpdate=RealDataValidation(app.MdataEngy,app.OldIds,app.SSPDsorted,app.NewIds);
    %app.MdataEngy=[];app.SSPDsorted=[];
    if ~isempty(MdataEngyUpdate)
        param=TableEstimatedParameters(1);
        
        CID=sum(~cellfun(@isempty,MdataEngyUpdate),2);
        EngBias=RealInspectorsEstimatedSigma;
        %%
        dt=1;
        A=[1 dt (dt^2)/2;0 1 dt;0 0 1];
        Q= param.^2*[(dt^5)/20 (dt^4)/8 (dt^3)/6;(dt^4)/8 (dt^3)/3 (dt^2)/2;(dt^3)/6 (dt^2)/2 dt];
        C=[1,0,0];
        % initial knowledge at time t=0
        init_x=[100;TableEstimatedParameters(5);TableEstimatedParameters(6)];
        % User Input
        M=10;init_V=zeros(3);DuplicatedIndex=[];RS=1;
        RandAcceptedData=randi([1 length(CID)],1,NumRealElements);
        CountVals=1;
        for i=1:length(CID)
            for ElementInd=1:CID(i)
                y=MdataEngyUpdate{i,ElementInd}(:,1);
                MaxDiff=max(abs(diff(y(1:end))));
                Insp=MdataEngyUpdate{i,ElementInd}(:,2)';
                if ~isempty(RegressionModel)
                    app.StructureElements=MdataEngyUpdate(i,1:CID(i));
                end
                if sum(isnan(y))<1 && length(find(~isnan(y)))>3 && ismember(Insp(end),EngBias(:,3)) && MaxDiff<=15 && 2*length(find(diff(y)>5))<length(diff(y))
                    if ismember(i,RandAcceptedData) && ElementInd==1 && MaxDiff<=15
                        app.BatchMode=get(app.PlotTimeSeries,'Value');%1
                    elseif MaxDiff>15
                        app.BatchMode=get(app.PlotTimeSeries,'Value');%1
                    else
                        app.BatchMode=0;
                    end
                    IndFinalObs=find(~isnan(y));
                    Yommited(CountVals,1)=y(IndFinalObs(end));
                    StructureInd=MdataEngyUpdate{i,ElementInd}(1,4);
                    yearly=MdataEngyUpdate{i,ElementInd}(1,3):MdataEngyUpdate{i,ElementInd}(end-1,3)+M;
                    DeltaYears(CountVals,1)=MdataEngyUpdate{i,ElementInd}(end,3)-MdataEngyUpdate{i,ElementInd}(end-1,3);
                    [~,~,iinsp]=intersect(MdataEngyUpdate{i,ElementInd}(:,3),yearly);
                    ObsYears=zeros(length(yearly),1);
                    ObsYears(iinsp)=1;
                    Diffy=nan(length(yearly),1);
                    if length(y)>length(iinsp)
                        [~, IndexValues] = unique(MdataEngyUpdate{i,ElementInd}(:,3), 'rows');
                        DuplicatedIndex=setdiff(1:size(MdataEngyUpdate{i,ElementInd}(:,3), 1), IndexValues);
                        y(DuplicatedIndex)=[];
                    end
                    Diffy(iinsp)=y;
                    yinsp=y;
                    y=Diffy';
                    % 0 is any inspector
                    % NaN is Unkown inspector (prediction)
                    Insp=MdataEngyUpdate{i,ElementInd}(:,2)';
                    ObsIndexes=find(~isnan(y));                                                 % Identify the indecies of observations (to take the average)
                    if length(ObsIndexes)==2
                        ObsIndexes=ObsIndexes(1:2);
                    elseif length(ObsIndexes)>2
                        ObsIndexes=ObsIndexes(1:3);
                    else
                        ObsIndexes=ObsIndexes(1);
                    end
                    AllAtt=[MdataEngyUpdate{i,ElementInd}(1,5:7) mean(y(ObsIndexes))];
                    if length(Insp)>length(find(~isnan(y)))
                        Insp(DuplicatedIndex)=[];
                    end
                    for k=1:length(Insp)
                        [~,ia,~] = intersect(EngBias(:,3),Insp(k));
                        if isempty(ia)
                            Re(k)=mean(EngBias(:,2)).^2;
                            InpecBiase(k)=0;%EngBias(ia,1);
                            Inspemmited=Insp(k);
                        else
                            Re(k)=EngBias(ia,2).^2;
                            InpecBiase(k)=0;%EngBias(ia,1);
                            Inspemmited=EngBias(ia,3);
                        end
                    end
                    Remmited(CountVals,1)=sqrt(Re(end));
                    InspOmmited(CountVals,1)=Inspemmited;
                    ReZero=zeros(1,length(y));
                    InpecBiaseZero=zeros(1,length(y));
                    InspectorLabel=nan(1,length(y));
                    ReZero(iinsp)=Re;
                    InpecBiaseZero(iinsp)=InpecBiase;
                    InspectorLabel(iinsp)=Insp;
                    Re=ReZero;
                    InpecBiase=InpecBiaseZero;
                    OptmInsp=zeros(1,length(y));
                    RU=1;InspBU=0;R=0;
                    % Excute
                    [loglik,~,~,~,~,~,~,~,~,~,~,~,~,x,s_Xsmooth,Vvalues]=KFsKF(app,y, A, C, Q, R, Re,...
                        init_x, init_V,InpecBiase,OptmInsp,RU,InspBU,ObsYears,yearly,...
                        InspectorLabel,StructureInd,ElementInd,TableEstimatedParameters(1,end),...
                        [],[],[],RegressionModel,AllAtt,TableEstimatedParameters,1,1);
                    
                    XInitialValuesStore(1:3,CountVals)=x(:,1);
                    XassociatedSE(1:2,CountVals)=[ElementInd;StructureInd];
                    IndFinalObs=find(~isnan(y))+1;
                    CondEstimate(CountVals,1)=x(1,IndFinalObs(end));
                    CondUncertainty(CountVals,1)=sqrt(Vvalues(1,1,IndFinalObs(end)));
                    
                    
                    StructuralAttValues=MdataEngyUpdate(i,1:CID(i));
                    CovariatesValues(CountVals,:)=[ElementInd StructuralAttValues{1,ElementInd}(1,5:end)];
                    
                    CountVals=CountVals+1;
                end
                clear y Re ReZero StructureInd yearly Insp RU InspBU R InpecBiase
            end
        end
        for i=1:length(InspOmmited)
            EstmInspector=find(EngBias(:,3)==InspOmmited(i,1));
            if ~isempty(EstmInspector)
                IncludedIndex(i)=i;
            end
        end
        CorrInd=find(IncludedIndex==0);
        IncludedIndex(CorrInd)=[];
        for i=1:length(IncludedIndex)
            EstimateCond(i,1)=(RevSpaceTransform(TableEstimatedParameters(1,end),...
                CondEstimate(IncludedIndex(i),1)));
            EstimateCondTr(i,1)=CondEstimate(IncludedIndex(i),1);
        end
        Markers={'o','*','s','x','p','h','+','d','.','<','>','^','v'};
        
        figure(1)
        HiddenObservation=Yommited(IncludedIndex,1);
        [~,HiddenObservationTr]=SpaceTransformation(TableEstimatedParameters(1,end),...
            Yommited(IncludedIndex,1),100.0001,25);
        DY=DeltaYears(IncludedIndex,1);
        UDY=unique(DY);
        for i=1:length(UDY)
            ScatterIndex=find(DY==UDY(i));
            scatter(EstimateCond(ScatterIndex),HiddenObservation(ScatterIndex),Markers{i});
            hold on
        end
        %scatter(EstimateCond,HiddenObservation,[],DeltaYears(IncludedIndex,1));
        xlabel('Model Estimate $\tilde{\mu}_{t|\mathtt{T}-1}$','Interpreter','latex')
        ylabel('Hidden Observation $\tilde{y}_{t=\mathtt{T}}$','Interpreter','latex')
        axis square
        grid on
        hold on
        plot([40,100],[40,100]);
        h=legend('1 Year','2 Years','3 Years','4 Years');
        set(h,'Interpreter','latex')
        set(h,'Position',[0.25 0.25 0.1 0.1])
        %colorbar
        hold off
        
        figure(2)
        Mu=EstimateCondTr-HiddenObservationTr';
        VarianceCond=CondUncertainty(IncludedIndex,1).^2+Remmited(IncludedIndex,1).^2;
        scatter(EstimateCond,HiddenObservation,(VarianceCond));%,[],sqrt(VarianceCond));
        for i=1:length(Mu)
            loglik(i)= gaussian_prob(Mu(i), zeros(1,length(Mu(i))), VarianceCond(i), 1);
        end
        app.LLEditField.Value=sum(loglik);
        xlabel('Model Estimate $\tilde{\mu}_{t|\mathtt{T}-1}$','Interpreter','latex')
        ylabel('Hidden Observation $\tilde{y}_{t=\mathtt{T}}$','Interpreter','latex')
        axis square
        grid on
        hold on
        plot([40,100],[40,100]);
        colorbar
        hold off
        
        figure(3)
        histogram(Mu./sqrt(VarianceCond),'normalization','pdf')
        hold on
        plot(-10:0.1:10,normpdf(-10:0.1:10,0,1))
        grid on
        title('Histogram: $[\mu_{t|\mathtt{T}-1}-y_{t=\mathtt{T}}]/\sigma_{t|\mathtt{T}-1}$','Interpreter','latex')
        hold off
        
%         figure(3)
%         DY=DeltaYears(IncludedIndex,1);
%         UDY=unique(DY);
%         for i=1:length(UDY)
%             ScatterIndex=find(DY==UDY(i));
%             scatter(EstimateCond(ScatterIndex),HiddenObservation(ScatterIndex),[],sqrt(VarianceCond(ScatterIndex)),Markers{i});
%             hold on
%         end
%         xlabel('Model Estimate $\tilde{\mu}_{t|\mathtt{T}-1}$','Interpreter','latex')
%         ylabel('Hidden Observation $\tilde{y}_{t=\mathtt{T}}$','Interpreter','latex')
%         axis square
%         grid on
%         box on
%         hold on
%         plot([40,100],[40,100]);
%         h=legend('1 Year','2 Years','3 Years','4 Years');
%         set(h,'Interpreter','latex')
%         set(h,'Position',[0.2 0.75 0.1 0.1])
%         colorbar
%         hold off
    else
        msgbox('There are no differences between the selected databases, analyses are not possible')
    end
        
end