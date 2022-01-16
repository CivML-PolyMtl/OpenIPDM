function [MdataEngy]=OrgnizeData(MdataEngy,CurrentInspectorID,...
    CurrentInspectorParam,EngBiasData,DataTransformedSpace,Ncurve,...
    OptBoundsData,AW,GPUCompute,IncludeStructuralAtt,NumStructuralAtt,TestSet)
rng(6);
R=0;IndexIt=[];jj=1;
CID=sum(~cellfun(@isempty,MdataEngy),2);
TrainSize=sum(CID);
TestingData=round(TestSet*TrainSize);  % 1200 for Tablier|3000 poutre Mur de front: 1500 % size of independent test set
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
                if max(abs(diff(y)))>15
                    IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                    y=nan(5,1);
                elseif length(diff(y))>2 && 2*length(find(diff(y)>5))>=length(diff(y))
                    IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                    y=nan(5,1);
                elseif length(diff(y))==2 && length(find(diff(y)>5))>=length(diff(y))
                    IndexIt=[IndexIt;MdataEngy{i,j}(1,4) j];
                    y=nan(5,1);
                elseif MdataEngy{i,j}(1,3)<2000 && DataTransformedSpace==0
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
                    ElementInd=j;
                end
            end
        else
            y=nan(5,1);
        end
        if sum(~isnan(y))>=3 
            if DataTransformedSpace==0
                [~,y]=SpaceTransformation(Ncurve,y,100,25);
            end
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
        end
        jj=jj+1;
    end
end

% check if the first obs is nan
NaNTimeSeries=isnan(YS(1,:,2));
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

MinInspCount=1;
while MinInspCount
    DataAllInspectors=AllInspectors;
    j=0;
    AllTestInd=[];
    while j<TestingData
        AllInspc=zeros(size(DataAllInspectors{1,1}));
        for i=1:length(DataAllInspectors)
            InspectorsCount(i)=sum(sum(DataAllInspectors{1,i}));
            AllInspc=AllInspc+DataAllInspectors{1,i}.*(1-1/InspectorsCount(i));
        end
        AllInspc((AllInspc==0))=nan;
        TSscore=min(AllInspc,[],2);
        TSscore(isnan(TSscore))=0;
        if sum(TSscore)~=0
            TestStrucIndex = randsample(1:length(AllInspc), 1, true, TSscore);
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
        end
    end
    if TestingData==0
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
    if MinCountVal~=0
        MinInspCount=0;
        break;
    end
end
TestStrucIndex=AllTestInd;
SetSize=round(length(TestStrucIndex)/3);
IndValidModel=0*SetSize+1:2*SetSize;
IndTestModel=2*SetSize+1:length(TestStrucIndex);
TestInd=TestStrucIndex(IndTestModel);
ValidInd=TestStrucIndex(IndValidModel);
TrainInd=1:size(AllInspectors{1,1},1);
if TestingData~=0
    TrainInd(TestStrucIndex)=[];
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

if TestingData==0
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
MdataEngy.AllInspectors=cellfun(@(x) x(TrainInd,:),...
    AllInspectors,'UniformOutput',false);
MdataEngy.init_x=init_x(1,TrainInd);
if IncludeStructuralAtt
    MdataEngy.StrucAtt=StrucAtt(1,TrainInd,:);
end

end