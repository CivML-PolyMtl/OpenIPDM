function InspectorStrucIndex=FindInspStrucIndex(MdataEngy,UID)
InspectorStrucIndex={};
for j=1:length(UID)
    InspectStrucIndex = cellfun(@(x) x==UID(j), MdataEngy, 'UniformOutput', 0);
    for i=1:length(InspectStrucIndex)
        SingleRow=InspectStrucIndex(i,:);
        V1=sum(cell2mat(cellfun(@sum,SingleRow,'un',0)));
        if V1>0
            InspectorStrucIndex{i,j}=i;
        end
    end
end
end