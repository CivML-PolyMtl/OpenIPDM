MdataEngy=SSPDsorted;
CID=sum(~cellfun(@isempty,MdataEngy),2);
fx=zeros(length(CID),max(CID));
Inspectors=[];
for i=1:length(CID)
    for j=1:CID(i)
        UniqueInspectors=unique(MdataEngy{i,j}(:,2));
        Inspectors = union(UniqueInspectors,Inspectors);
    end
end

InspObs={};InspStruc={};ObsStore=[];StrucStore=[];
for k=1:length(Inspectors)
    for i=1:length(CID)
        for j=1:CID(i)
            InspecIndex=find(MdataEngy{i,j}(:,2)==Inspectors(k));
            if ~isempty(InspecIndex)
                ObsStore=[ObsStore MdataEngy{i,j}(InspecIndex,1)'];
                StrucStore=[StrucStore MdataEngy{i,j}(InspecIndex,4)'];
            end
        end
    end
    InspObs{k}=ObsStore;
    InspStruc{k}=StrucStore;
    ObsStore=[];StrucStore=[];
end
Vc = cellfun(@(x) unique(x(:)), InspStruc, 'UniformOutput',false);
NumStruc=cellfun('size',Vc,1);
ObsMean = cellfun(@(x) nanmean(x(:)), InspObs, 'UniformOutput',false);
ObsStd = cellfun(@(x) nanstd(x(:)), InspObs, 'UniformOutput',false);
NumObs=cellfun('size',InspObs,2);
Inspectors=[Inspectors NumStruc' NumObs' cell2mat(ObsMean)' cell2mat(ObsStd)'];