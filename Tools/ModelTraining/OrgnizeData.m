function [MdataEngy, Success]=OrgnizeData(MdataEngy,CurrentInspectorID,...
    CurrentInspectorParam,EngBiasData,AW,GPUCompute,IncludeStructuralAtt,...
    NumStructuralAtt,TestSet,IncludeInterventions, syn_flag)
rng(6);
R=0;IndexIt=[];jj=1;
CID=sum(~cellfun(@isempty,MdataEngy),2);
TrainSize=sum(CID);
TestingData=round(TestSet*TrainSize);  % 1200 for Tablier|3000 poutre Mur de front: 1500 % size of independent test set
Success = 1;
if GPUCompute==2
    YS=NaN(1,TrainSize,AW,'gpuArray');
    RS=zeros(1,TrainSize,1,'gpuArray');
    ReS=zeros(1,TrainSize,AW,'gpuArray');
    InpecBiaseS=zeros(1,TrainSize,AW,'gpuArray');
    CurrentInspectorS=zeros(1,TrainSize,AW,'gpuArray');
    RUS=zeros(1,TrainSize,AW,'gpuArray');
    InspBUS=zeros(1,TrainSize,AW,'gpuArray');
    ObsYearsS=zeros(1,TrainSize,AW,'gpuArray');
    yearlyS=zeros(1,TrainSize,AW,'gpuArray');
    InspectorLabelS=NaN(1,TrainSize,AW,'gpuArray');
    StructureIndS=zeros(1,TrainSize,1,'gpuArray');
    AllInspectors=cell(1,size(EngBiasData,1));
    ElementIndS=zeros(1,TrainSize,1,'gpuArray');
    if IncludeStructuralAtt
        StrucAtt=zeros(1,TrainSize,NumStructuralAtt,'gpuArray');
    end
    if IncludeInterventions
        InterventionYear=zeros(1,TrainSize,1,'gpuArray');
        InterventionType=zeros(1,TrainSize,1,'gpuArray');
        InterventionJump=zeros(1,TrainSize,1,'gpuArray');
        InterventionVector=zeros(1,TrainSize,AW,'gpuArray');
    end    
else
    YS=NaN(1,TrainSize,AW);
    RS=zeros(1,TrainSize,1);
    ReS=zeros(1,TrainSize,AW);
    InpecBiaseS=zeros(1,TrainSize,AW);
    CurrentInspectorS=zeros(1,TrainSize,AW);
    RUS=zeros(1,TrainSize,AW);
    InspBUS=zeros(1,TrainSize,AW);
    ObsYearsS=zeros(1,TrainSize,AW);
    yearlyS=zeros(1,TrainSize,AW);
    InspectorLabelS=NaN(1,TrainSize,AW);
    StructureIndS=zeros(1,TrainSize,1);
    AllInspectors=cell(1,size(EngBiasData,1));
    ElementIndS=zeros(1,TrainSize,1);
    if IncludeStructuralAtt
        StrucAtt=zeros(1,TrainSize,NumStructuralAtt);
    end
    if IncludeInterventions
        InterventionYear=zeros(1,TrainSize,AW);
        InterventionType=zeros(1,TrainSize,1);
        InterventionJump=zeros(1,TrainSize,1);
        InterventionVector=zeros(1,TrainSize,AW);
    end
end


for i=1:length(CID)
    for j=1:CID(i)
        Re=[];InpecBiase=[];
        % stretch the y vector to 1 year basis
        while length(unique(MdataEngy{i,j}(:,3)))~=length(MdataEngy{i,j}(:,3))
            [~, ii] = unique(MdataEngy{i,j}(:,3), 'first');
            keep    = (diff([ii(:); length(MdataEngy{i,j}(:,3)) + 1]) > 1);
            Repindex  = ii(keep);
            if ~isempty(Repindex)
                RepInd=Repindex+1;
                MdataEngy{i,j}(RepInd,:)=[];
            end
        end
        NanDatabase=find(isnan(MdataEngy{i,j}(:,1)'));
        MdataEngy{i,j}(NanDatabase,:)=[];
        if ~isempty(MdataEngy{i,j})
            y=MdataEngy{i,j}(:,1)'; %y (database)
            if length(y)<3
                y=nan(5,1);
                IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
            else
                if max(abs(diff(y)))>15 && ~IncludeInterventions
                  IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                  y=nan(5,1);
                elseif length(diff(y))>2 && 3*length(find(diff(y)>10))>=length(diff(y)) && ~IncludeInterventions
                  IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                  y=nan(5,1);
                elseif length(diff(y))>2 && 2*length(find(diff(y)>5))>=length(diff(y)) && ~IncludeInterventions
                  IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                  y=nan(5,1);
                elseif length(diff(y))>2 && all(diff(y(~isnan(y)))>5) && ~IncludeInterventions
                  IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                  y=nan(5,1);
                elseif length(diff(y))==2 && length(find(diff(y)>5))>=length(diff(y)) && ~IncludeInterventions
                  IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                  y=nan(5,1);
                elseif MdataEngy{i,j}(1,3)<2000 && ~syn_flag
                    IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                    y=nan(5,1);
                else
                    % generate time sequence
                    yearly=MdataEngy{i,j}(1,3):MdataEngy{i,j}(end,3);
                    
                    % find the years where that have observations
                    [~,~,iinsp]=intersect(MdataEngy{i,j}(:,3),yearly);
                    ObsYears=zeros(length(yearly),1);
                    ObsYears(iinsp)=1; 
                    % stretch the observation vector over the time series
                    Diffy=nan(length(yearly),1);
                    Diffy(iinsp)=y;
                    y=Diffy';
                    % find inspectors for the time series
                    Insp=MdataEngy{i,j}(:,2)';
                    % 0 is any inspector
                    % NaN is Unkown inspector
                    % 1 is current inspector (to be optimized)
                    CurrentInspector=zeros(1,length(y));
                    CurrentInspector(find(isnan(y)))=NaN;
                    AInsp=zeros(1,length(y));
                    AInsp(iinsp)=Insp;
                    [~,~,ib] = intersect(CurrentInspectorID,AInsp);
                    CurrentInspector(ib)=1;
                    % the (to be optimized) inspector biase
                    RUv=zeros(1,length(y)); 
                    RU=(CurrentInspectorParam(1)).^2;
                    RUv(ib)=RU;
                    RU=RUv;
                    InspBU=0;%CurrentInspectorParam(1);
                    % fill the inspector biase and inspector variance for the
                    % time series
                    InspIndexV=find(~isnan(CurrentInspector));
                    for k=1:length(Insp)
                        [~,ia,~] = intersect(EngBiasData(:,1),Insp(k));
                        if isempty(ia)
                            EngBiasData=[EngBiasData;[Insp(k) EngBiasData(end,2:3)]];
                            AllInspectors=[AllInspectors cell(1,1)];
                            [~,ia,~] = intersect(EngBiasData(:,1),Insp(k));
                        end
                        Re(k)=(EngBiasData(ia,end)).^2;
                        InpecBiase(k)=0;%EngBiasData(ia,1);
                    end
                    ReZero=zeros(1,length(y));
                    InpecBiaseZero=zeros(1,length(y));
                    InspectorLabel=nan(1,length(y));
                    ReZero(iinsp)=Re;
                    InpecBiaseZero(iinsp)=InpecBiase;
                    InspectorLabel(iinsp)=Insp;
                    Re=ReZero;
                    InpecBiase=InpecBiaseZero;
                    StructureInd=MdataEngy{i,j}(1,4);
                    if IncludeStructuralAtt
                        StructuralAttributeValue=MdataEngy{i,j}(1,5:5+NumStructuralAtt-1);
                    end
                    if IncludeInterventions
                        AllInterType=unique(MdataEngy{i,j}(1:end,end));
                        AllInterTypeIndex=find(AllInterType>0);
                        AllInterYear=unique(MdataEngy{i,j}(1:end,end-1));
                        AllInterYearIndex=find(AllInterYear>0);
                        if ~isempty(AllInterTypeIndex) && ~isempty(AllInterYearIndex)
                            InterType=AllInterType(AllInterTypeIndex(end));
                            InterType=fix(InterType./10.^fix(log10(InterType)));
                            InterYear=AllInterYear(AllInterYearIndex(end));
                            if InterType==4
                                y=nan(5,1);
                            end
                            
                            if InterYear~=0
                                IndexLast=find(InterYear<=MdataEngy{i,j}(:,3));
                                InterventionVectorElm=zeros(1,length(yearly)+1);
                                IntIndYear=find(yearly==InterYear);
                                if ~isempty(IntIndYear) && IntIndYear>2 %&& isempty(find(InterYear==MdataEngy{i,j}(:,3)))
%                                    if DataTransformedSpace==0
                                       InterventionVectorElm(IntIndYear+1)=1;
%                                    else
%                                        InterventionVectorElm(IntIndYear-1)=1;
%                                    end
%                                 elseif ~isempty(IntIndYear) && IntIndYear>2 && ~isempty(find(InterYear==MdataEngy{i,j}(:,3)))
%                                     InterventionVectorElm(IntIndYear+1)=1;
                                else
                                    IndexLast=[];
                                    y=nan(5,1);
                                end
                            else
                                IndexLast=[];
                                y=nan(5,1);
                            end
                            if ~isempty(IndexLast)
                                if IndexLast(1)>1
                                    InterJump=MdataEngy{i,j}(IndexLast(1),1)-MdataEngy{i,j}(IndexLast(1)-1,1);
                                else
                                    InterType=0;
                                    InterYear=0;
                                    InterJump=0;
                                    y=nan(5,1);
                                end
                                if InterType>2 && MdataEngy{i,j}(IndexLast(1)-1,1)>90
                                    InterType=0;
                                    InterYear=0;
                                    InterJump=0;
                                    y=nan(5,1);
                                end
                            else
                                InterType=0;
                                InterYear=0;
                                InterJump=0;
                                y=nan(5,1);
                            end
                        else
                            InterType=0;
                            InterYear=0;
                            InterJump=0;
                            y=nan(5,1);
                        end
                        if InterType>1 && InterJump<=0.1
                            InterType=0;
                            InterYear=0;
                            InterJump=0;
                            y=nan(5,1);
                        end
                        
                        if InterJump>50
                            InterType=0;
                            InterYear=0;
                            InterJump=0;
                            y=nan(5,1);
                        end
                    end
                    ElementInd=j;
                end
            end
        else
            y=nan(5,1);
        end
        if sum(~isnan(y))>=3 
            YS(1,jj,2:length(y)+1)=y;
            RS(1,jj,:)=R;
            ReS(1,jj,2:length(Re)+1)=Re;
            InpecBiaseS(1,jj,2:length(InpecBiase)+1)=InpecBiase;
            CurrentInspectorS(1,jj,2:length(CurrentInspector)+1)=CurrentInspector;
            RUS(1,jj,2:length(RU)+1)=RU;
            InspBUS(1,jj,:)=InspBU;
            ObsYearsS(1,jj,2:length(ObsYears)+1)=ObsYears;
            yearlyS(1,jj,1:length(yearly))=yearly;
            InspectorLabelS(1,jj,2:length(InspectorLabel)+1)=InspectorLabel;
            NaNIndex=find(~isnan(InspectorLabel));
            UILabel=unique(InspectorLabel(NaNIndex));
            for II=1:length(UILabel)
                ObsInspIdx=find(InspectorLabel==UILabel(II));
                InspectorIndex=find(EngBiasData(:,1)==UILabel(II));
                if isempty(AllInspectors{InspectorIndex})
                    AllInspectors{InspectorIndex}=zeros(TrainSize,AW);
                end
                AllInspectors{InspectorIndex}(jj,ObsInspIdx+1)=1;
            end
            StructureIndS(1,jj,:)=StructureInd;
            ElementIndS(1,jj,:)=ElementInd;
            if IncludeStructuralAtt
                StrucAtt(1,jj,:)=StructuralAttributeValue;
            end
            if IncludeInterventions
                InterventionYear(1,jj,1)=InterYear;
                InterventionType(1,jj,1)=InterType;
                InterventionJump(1,jj,1)=InterJump;
                InterventionVector(1,jj,1:length(InterventionVectorElm))=InterventionVectorElm;
            end
        end
        jj=jj+1;
    end
end

% check if the first obs is nan
NaNTimeSeries=isnan(YS(1,:,2));
if sum(NaNTimeSeries)~= length(YS(1,:,2))
    
    if ~isempty(NaNTimeSeries)
        YS(:,NaNTimeSeries,:)=[];
        RS(:,NaNTimeSeries)=[];
        ReS(:,NaNTimeSeries,:)=[];
        InpecBiaseS(:,NaNTimeSeries,:)=[];
        CurrentInspectorS(:,NaNTimeSeries,:)=[];
        RUS(:,NaNTimeSeries,:)=[];
        InspBUS(:,NaNTimeSeries,:)=[];
        ObsYearsS(:,NaNTimeSeries,:)=[];
        yearlyS(:,NaNTimeSeries,:)=[];
        InspectorLabelS(:,NaNTimeSeries,:)=[];
        StructureIndS(:,NaNTimeSeries,:)=[];
        ElementIndS(:,NaNTimeSeries)=[];
        if IncludeStructuralAtt
            StrucAtt(:,NaNTimeSeries,:)=[];
        end
        if IncludeInterventions
            InterventionYear(:,NaNTimeSeries)=[];
            InterventionType(:,NaNTimeSeries)=[];
            InterventionJump(:,NaNTimeSeries)=[];
            InterventionVector(:,NaNTimeSeries,:)=[];
        end
        for i=1:length(AllInspectors)
            if ~isempty(AllInspectors{i})
                AllInspectors{i}(NaNTimeSeries,:)=[];
            end
        end
    end

    % compute the average of the first 2:3 observations
    init_x=zeros(1,length(YS(1,:,2)));
    for i=1:length(YS(1,:,2))
        IndexInit=find(~isnan(gather(YS(1,i,:))));
        if length(IndexInit)==2
            IndexInit=IndexInit(1:2);
        elseif length(IndexInit)>2
            IndexInit=IndexInit(1:3);
        else
            IndexInit=IndexInit(1);
        end
        init_x(1,i)=mean(gather(YS(1,i,IndexInit)));
    end

    TimedOut = 0;
    MinInspCount=1;
    while MinInspCount
        DataAllInspectors=AllInspectors;
        j=0;
        AllTestInd=[];
        while j<TestingData
            AllInspc=zeros(size(DataAllInspectors{1,1}));
            for i=1:length(DataAllInspectors)
                InspectorsCount(i)=sum(sum(DataAllInspectors{1,i}));
                % any error in this section could be connected to a
                % missmatch in the filtering criteria of the training data
                % between the function "FIlteredInspectors" and THIS function
                AllInspc=AllInspc+DataAllInspectors{1,i}.*(1-1/InspectorsCount(i));
            end
            AllInspc((AllInspc==0))=nan;
            TSscore=min(AllInspc,[],2);
            TSscore(isnan(TSscore))=0;
            if sum(TSscore)~=0
                TestStrucIndex = randsample(1:size(AllInspc,1), 1, true, TSscore);
                %%
                TestIndStruc=find(StructureIndS==StructureIndS(TestStrucIndex));
                AllTestInd=[AllTestInd;TestIndStruc'];
                for i=1:length(DataAllInspectors)
                    DataAllInspectors{1,i}(AllTestInd,:)=0;
                end
                j=length(AllTestInd);
            else
                j=0;
                AllTestInd=[];
                DataAllInspectors=AllInspectors;
                if TimedOut > 50000
                    TestingData = 0;
                    break;
                end
                TimedOut = TimedOut + 1;
            end
        end
        if TestingData==0 && ~IncludeInterventions
            AllInspc=zeros(size(DataAllInspectors{1,1}));
            for i=1:length(DataAllInspectors)
                InspectorsCount(i)=sum(sum(DataAllInspectors{1,i}));
                AllInspc=AllInspc+DataAllInspectors{1,i}.*(1-1/InspectorsCount(i));
            end
            AllInspc((AllInspc==0))=nan;
            BinInspec=~isnan(AllInspc);
            NumObsTS=sum(BinInspec,2);
            NormNOT=NumObsTS./max(NumObsTS);
            LRmat=fliplr(BinInspec);
            [~,inmax]=max(LRmat,[],2);
            LastIndex=size(BinInspec,2)-inmax+1;
            CS=repmat(size(BinInspec,2),length(LastIndex)-1,1);
            CS=[0;cumsum(CS)];
            AllInspcTranspose=AllInspc';
            TSscore=AllInspcTranspose(CS+LastIndex);
            NormTS=TSscore./max(TSscore);
            NormTS(isnan(NormTS))=0;
            if sum(NormTS)~=0
                [~,SortedScore]=sort(NormTS+NormNOT,'descend');
                AllTestInd = SortedScore(1:round(2/3*size(BinInspec,1)));
            end
            break;
        end
        for i=1:length(DataAllInspectors)
            InspectorsCount(i)=sum(sum(DataAllInspectors{1,i}));
        end
        MinCountVal=min(InspectorsCount);
        if MinCountVal~=0 || IncludeInterventions
            MinInspCount=0;
            break;
        end
        if length(YS(1,:,2))<100
            TestingData = 0;
        end
        if TestingData > round(TestSet*length(YS(1,:,2)))
            TestingData = round(TestSet*length(YS(1,:,2)));
        end
    end

    TestStrucIndex=AllTestInd;
    SetSize=round(length(TestStrucIndex)/3);
    IndValidModel=0*SetSize+1:2*SetSize;
    if isempty(IndValidModel) && ~isempty(TestStrucIndex)
        IndValidModel = 1;
    end
    IndTestModel=2*SetSize+1:length(TestStrucIndex);
    TestInd=TestStrucIndex(IndTestModel);
    ValidInd=TestStrucIndex(IndValidModel);
    TrainInd=1:size(AllInspectors{1,1},1);
    if TestingData~=0
        TrainInd(TestStrucIndex)=[];
    end

    if TestingData==0 && IncludeInterventions
        TrainInd = 1:length(YS(1,:,1));
    end
    MdataEngy=[];

    MdataEngy.ModelValid.YS=YS(1,ValidInd,:);
    MdataEngy.ModelValid.RS=RS(1,ValidInd);
    MdataEngy.ModelValid.ReS=ReS(1,ValidInd,:);
    MdataEngy.ModelValid.InpecBiaseS=InpecBiaseS(1,ValidInd,:);
    MdataEngy.ModelValid.CurrentInspectorS=CurrentInspectorS(1,ValidInd,:);
    MdataEngy.ModelValid.RUS=RUS(1,ValidInd,:);
    MdataEngy.ModelValid.InspBUS=InspBUS(1,ValidInd,:);
    MdataEngy.ModelValid.ObsYearsS=ObsYearsS(1,ValidInd,:);
    MdataEngy.ModelValid.yearlyS=yearlyS(1,ValidInd,:);
    MdataEngy.ModelValid.InspectorLabelS=InspectorLabelS(1,ValidInd,:);
    MdataEngy.ModelValid.StructureIndS=StructureIndS(1,ValidInd);
    MdataEngy.ModelValid.ElementIndS=ElementIndS(1,ValidInd);
    MdataEngy.ModelValid.AllInspectors=cellfun(@(x) x(ValidInd,:),...
        AllInspectors,'UniformOutput',false);
    MdataEngy.ModelValid.init_x=init_x(1,ValidInd);
    if IncludeStructuralAtt
        MdataEngy.ModelValid.StrucAtt=StrucAtt(1,ValidInd,:);
    end
    if IncludeInterventions && ~isempty(ValidInd)
        MdataEngy.ModelValid.InterventionYear=InterventionYear(1,ValidInd);
        MdataEngy.ModelValid.InterventionType=InterventionType(1,ValidInd);
        MdataEngy.ModelValid.InterventionJump=InterventionJump(1,ValidInd);
        MdataEngy.ModelValid.InterventionVector=InterventionVector(1,ValidInd,:);
    end

    MdataEngy.ModelTest.YS=YS(1,TestInd,:);
    MdataEngy.ModelTest.RS=RS(1,TestInd);
    MdataEngy.ModelTest.ReS=ReS(1,TestInd,:);
    MdataEngy.ModelTest.InpecBiaseS=InpecBiaseS(1,TestInd,:);
    MdataEngy.ModelTest.CurrentInspectorS=CurrentInspectorS(1,TestInd,:);
    MdataEngy.ModelTest.RUS=RUS(1,TestInd,:);
    MdataEngy.ModelTest.InspBUS=InspBUS(1,TestInd,:);
    MdataEngy.ModelTest.ObsYearsS=ObsYearsS(1,TestInd,:);
    MdataEngy.ModelTest.yearlyS=yearlyS(1,TestInd,:);
    MdataEngy.ModelTest.InspectorLabelS=InspectorLabelS(1,TestInd,:);
    MdataEngy.ModelTest.StructureIndS=StructureIndS(1,TestInd);
    MdataEngy.ModelTest.ElementIndS=ElementIndS(1,TestInd);
    MdataEngy.ModelTest.AllInspectors=cellfun(@(x) x(TestInd,:),...
        AllInspectors,'UniformOutput',false);
    MdataEngy.ModelTest.init_x=init_x(1,TestInd);
    if IncludeStructuralAtt
        MdataEngy.ModelTest.StrucAtt=StrucAtt(1,TestInd,:);
    end
    if IncludeInterventions && ~isempty(TestInd)
        MdataEngy.ModelTest.InterventionYear=InterventionYear(1,TestInd);
        MdataEngy.ModelTest.InterventionType=InterventionType(1,TestInd);
        MdataEngy.ModelTest.InterventionJump=InterventionJump(1,TestInd);
        MdataEngy.ModelTest.InterventionVector=InterventionVector(1,TestInd,:);
    end

    if TestingData==0 && ~IncludeInterventions
        LLI=LastIndex(TestStrucIndex);
        for i=1:length(TestStrucIndex)
            YS(1,TestStrucIndex(i),LLI(i))=nan;
            RS(1,TestStrucIndex(i),LLI(i))=nan;
            ReS(1,TestStrucIndex(i),LLI(i))=nan;
            InspBUS(1,TestStrucIndex(i),LLI(i))=nan;
            InpecBiaseS(1,TestStrucIndex(i),LLI(i))=nan;
            ObsYearsS(1,TestStrucIndex(i),LLI(i))=0;
            StoreInspLabel=InspectorLabelS(1,TestStrucIndex(i),LLI(i));
            InspectorIndexTest=find(EngBiasData(:,1)==StoreInspLabel);
            AllInspectors{InspectorIndexTest}(TestStrucIndex(i),LLI(i))=0;
        end
    end

    MdataEngy.RemovedData=IndexIt;
    MdataEngy.YS=YS(1,TrainInd,:);
    MdataEngy.RS=RS(1,TrainInd);
    MdataEngy.ReS=ReS(1,TrainInd,:);
    MdataEngy.InpecBiaseS=InpecBiaseS(1,TrainInd,:);
    MdataEngy.CurrentInspectorS=CurrentInspectorS(1,TrainInd,:);
    MdataEngy.RUS=RUS(1,TrainInd,:);
    MdataEngy.InspBUS=InspBUS(1,TrainInd,:);
    MdataEngy.ObsYearsS=ObsYearsS(1,TrainInd,:);
    MdataEngy.yearlyS=yearlyS(1,TrainInd,:);
    MdataEngy.InspectorLabelS=InspectorLabelS(1,TrainInd,:);
    MdataEngy.StructureIndS=StructureIndS(1,TrainInd);
    MdataEngy.ElementIndS=ElementIndS(1,TrainInd);
    if TestingData==0 && ~IncludeInterventions
        MdataEngy.AllInspectors=cellfun(@(x) x(TrainInd,:),...
            AllInspectors,'UniformOutput',false);
    elseif TestingData==0 && IncludeInterventions
        MdataEngy.AllInspectors = AllInspectors;
    else
        MdataEngy.AllInspectors=cellfun(@(x) x(TrainInd,:),...
            AllInspectors,'UniformOutput',false);
    end
    MdataEngy.init_x=init_x(1,TrainInd);
    if IncludeStructuralAtt
        MdataEngy.StrucAtt=StrucAtt(1,TrainInd,:);
        MdataEngy.StrucAttCategories=unique(StrucAtt(1,:,1));
    end
    if IncludeInterventions
        MdataEngy.InterventionYear=InterventionYear(1,TrainInd);
        MdataEngy.InterventionType=InterventionType(1,TrainInd);
        MdataEngy.InterventionJump=InterventionJump(1,TrainInd);
        MdataEngy.InterventionVector=InterventionVector(1,TrainInd,:);
    end
else
     MdataEngy = [];
     Success = 0;
%     MdataEngy.RemovedData=NaNTimeSeries;
%     MdataEngy.YS=YS;
%     MdataEngy.RS=RS;
%     MdataEngy.ReS=ReS;
%     MdataEngy.InpecBiaseS=InpecBiaseS;
%     MdataEngy.CurrentInspectorS=CurrentInspectorS;
%     MdataEngy.RUS=RUS;
%     MdataEngy.InspBUS=InspBUS;
%     MdataEngy.ObsYearsS=ObsYearsS;
%     MdataEngy.yearlyS=yearlyS;
%     MdataEngy.InspectorLabelS=InspectorLabelS;
%     MdataEngy.StructureIndS=StructureIndS;
%     MdataEngy.ElementIndS=ElementIndS;
%     MdataEngy.AllInspectors=cellfun(@(x) x(:,:),...
%         AllInspectors,'UniformOutput',false);
%     MdataEngy.init_x=init_x;
%     if IncludeStructuralAtt
%         MdataEngy.StrucAtt=StrucAtt;
%     end
%     if IncludeInterventions
%         MdataEngy.InterventionYear=InterventionYear;
%         MdataEngy.InterventionType=InterventionType;
%         MdataEngy.InterventionJump=InterventionJump;
%         MdataEngy.InterventionVector=InterventionVector;
%     end

end