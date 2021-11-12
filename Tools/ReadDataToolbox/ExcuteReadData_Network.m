StatusCond=find(strcmp(app.SelectStructureActiveStatus.Items,app.SelectStructureActiveStatus.Value));
StructuresDatabaseStore.SelectedFormats(strcmp(StructuresDatabaseStore.SelectedFormats,'%f'))={'%q'};
StructuresData=readall(StructuresDatabaseStore);
if ~isempty(StatusCond)
    StatusColsVal=get(app.SelectStructureActiveStatus,'value');
    StatusColsInd=find(strcmp(app.SelectStructureActiveCol.Items,app.SelectStructureActiveCol.Value))-1;
    StatusRows=find(strcmp(StatusColsVal,table2array(StructuresData(:,StatusColsInd))));
    StrucVars=app.SelectStructureActiveCol.Items;
    StrucVars(1)=[];
    StructuresData=StructuresData(StatusRows,StrucVars);
end
StructuresData=table2array(StructuresData);
DJMAIndex=find(strcmp(app.SelectTrafficCol.Value,app.SelectTrafficCol.Items))-1;
TrucksIndex=find(strcmp(app.SelectTruckCol.Value,app.SelectTruckCol.Items))-1;
BLengthIndex=find(strcmp(app.SelectLengthCol.Value,app.SelectLengthCol.Items))-1;
NetworkData=StructuresData(:,1);
NetworkData_ElementDetails=StructuresData(:,1);
NetworkData(:,2)=StructuresData(:,DJMAIndex);
% NetworkData(:,3)=StructuresData(:,TrucksIndex);
% NetworkData(:,4)=StructuresData(:,BLengthIndex);
SEDatabaseStore.SelectedVariableNames=SelectedVars;
SEDatabaseStore.SelectedFormats(~strcmp(SEDatabaseStore.SelectedFormats,'%q'))={'%q'};
SEData=readall(SEDatabaseStore);
AllIDs = unique(SEData(:,1));
NewID={};
NewID(:,2)=table2cell(AllIDs);
NewID(:,1)=num2cell(1:height(AllIDs));
ElementDetailsDataStore.SelectedFormats(~strcmp(ElementDetailsDataStore.SelectedFormats,'%q'))={'%q'};
EDQ=readall(ElementDetailsDataStore);
FilterCond=find(strcmp(app.SelectStructuralElement.Items,app.SelectStructuralElement.Value));
for ElementType=1:length(PrimaryElements)
    NoIDIntRows=[];
    if FilterCond~=1
        ColsVal=PrimaryElements{ElementType};
        ColsInd=find(strcmp(app.SelectStructuralElement.Items,app.SelectStructuralElement.Value))-1;
        Rows=find(strcmp(ColsVal,table2array(SEData(:,ColsInd))));
        Vars=SelectedVars;
        SEDataEC=SEData(Rows,Vars);
        InterventionsRows=find(strcmp(ColsVal,table2array(InterventionSEData(:,ColsInd))));
        InterventionSEDataEL=InterventionSEData(InterventionsRows,[1 2 3 5]);
    end
    InterventionCond=find(strcmp(app.SelectStructureStatus.Items,app.SelectStructureStatus.Value));
    if InterventionCond~=1
        AllRows=[];
        InterAllRows=[];
        UISEData=unique(InterventionSEDataEL);
        ISID=str2double(table2cell(SEDataEC(:,end)));
        ISID_ind=str2double(table2cell(SEDataEC(:,1)));
        for i=1:size(table2cell(UISEData),1)
            SRows = find(strcmp(table2cell(UISEData(i,end)),table2cell(SEDataEC(:,end))));
            if isempty(SRows)
                SRows=find(ISID==str2double(cell2mat(table2cell(UISEData(i,end)))));
            end
            if isempty(SRows)
                SRows = find(strcmp(table2cell(UISEData(i,1)),table2cell(SEDataEC(:,1))));
                if isempty(SRows)
                    SRows=find(ISID_ind==str2double(cell2mat(table2cell(UISEData(i,1)))));
                end
                if ~isempty(SRows)
                    NoIDIntRows=[NoIDIntRows;table2cell(UISEData(i,1))];
                end
            end
            AllRows=[AllRows;SRows];
            if ~isempty(SRows)
                InterAllRows=[InterAllRows;i];
            end
        end
    end
    if InterventionCond==1
        SEDataEC(:,end)=[];
    end
    if InterventionCond==2
        NoInterRows=setdiff(1:size(SEDataEC,1),AllRows)';
        SEDataEC=SEDataEC(NoInterRows,:);
        SEDataEC(:,end)=[];
    elseif InterventionCond==3
        %     InterRows=intersect(1:size(SEData,1),AllRows)';
        SEDataEC=SEDataEC(AllRows,:);
        UISEData=UISEData(InterAllRows,:);
        SEDataEC(:,end)=[];
    end
    if ~isempty(EDQ) % Element details dataset: Struc ID, No Elem, No trav, Elem, QTE
        ElemRow=find(strcmp(ColsVal,table2array(EDQ(:,3))));
        EDQ_sub=table2array(EDQ(ElemRow,[1 2 4 5]));
        StructuralElemID=3;
        StructuralElemNumSpans=2;
        for i=1:size(NetworkData_ElementDetails,1)
            %% Unique elements
            ElemID=find(strcmp(NetworkData_ElementDetails(i,1),EDQ_sub(:,1)));
            if ~isempty(ElemID)
               QuantityData= str2double(EDQ_sub(ElemID,end));
               NoData=find(isnan(QuantityData));
               if ~isempty(NoData)
                   QuantityData(NoData)=nanmean(str2double(EDQ_sub(ElemID,end)));
                   NoDataAgain=find(isnan(QuantityData));
                   if ~isempty(NoDataAgain)
                       QuantityData(NoData)=nanmean(str2double(EDQ_sub(:,end)));
                       NetworkData_ElementDetails{i,ElementType+1}=QuantityData;
                   else
                        NetworkData_ElementDetails{i,ElementType+1}=QuantityData;
                   end
               else
                    NetworkData_ElementDetails{i,ElementType+1}=str2double(EDQ_sub(ElemID,end));
               end
            end
        end
    end
    
    InspectorsDatabase=readall(InspectorsDatabaseStore);
    InspectorsDatabase=table2array(InspectorsDatabase);
    SEDataEC=table2array(SEDataEC);
    EmptyEng=[];
    for i=1:size(SEDataEC,1)
        EngId=find(strcmp(SEDataEC(i,end),InspectorsDatabase(:,1)));
        if ~isempty(EngId)
            SEDataEC(i,end)=num2cell(InspectorsDatabase(EngId,end));
        else
            EmptyEng=[EmptyEng;i];
        end
    end
    
    SEDataEC(EmptyEng,:)=[];
    EndIndex=size(SEDataEC,2);
    for i=1:length(StructuresData)
        Ind=find(strcmp(StructuresData(i,1),SEDataEC(:,1)));
        for j=1:size(StructuresData,2)-1
            if ~isempty(Ind)
                SEDataEC(Ind,EndIndex+j)=StructuresData(i,j+1); % Const Date
            end
        end
    end
    
    ConstructionDateCol=find(strcmp(app.ConstYearCol.Items,app.ConstYearCol.Value))+EndIndex-1;
    InspectionYearCol=find(strcmp(app.InspecYearCol.Items,app.InspecYearCol.Value))-1;
    StructuralElemID=find(strcmp(app.SelectSElementID.Items,app.SelectSElementID.Value))-1;
    StructuralElemNumSpans=find(strcmp(app.SelectSNSpans.Items,app.SelectSNSpans.Value))-1;
    MaterialData=find(strcmp(app.SelectMaterial.Items,app.SelectMaterial.Value))-1;
    STAtt1=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt1.Value))+EndIndex-2;
    STAtt2=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt2.Value))+EndIndex-2;
    STAtt3=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt3.Value))+EndIndex-2;
    STAtt4=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt4.Value))+EndIndex-2;
    STAtt5=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt5.Value))+EndIndex-2;
    STAtt6=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt6.Value))+EndIndex-2;
    STAtt7=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt7.Value))+EndIndex-2;
    STAtt8=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt8.Value))+EndIndex-2;
    STAtt9=find(strcmp(app.SelectSAtt1.Items,app.SelectSAtt9.Value))+EndIndex-2;
    
    if MaterialData~=0
        MaterialCategories=unique(SEDataEC(:,MaterialData));
        if ~isempty(MaterialCategories)
            if ~isnumeric(MaterialCategories(1))
                for i=1:length(MaterialCategories)
                    SEDataEC(find(strcmp(MaterialCategories(i),SEDataEC(:,MaterialData))),MaterialData)={sprintf('%d',i)};
                end
            end
        end
    end
    
    for i=1:size(SEDataEC,1)
        DateVec=datestr(datenum(SEDataEC(i,InspectionYearCol)));
        SEDataEC(i,InspectionYearCol)={DateVec(end-3:end)};
    end
    
    EndIndex2=size(SEDataEC,2);
    for i=1:length(StructuresData(:,1))
        SRows = find(strcmp(StructuresData(i,1),SEDataEC(:,1)));
        SEDsorted{i,1}=SEDataEC(SRows,:);
        %% Add time index
        UniqueDate=unique(SEDsorted{i,1}(:,InspectionYearCol));
        if ~isempty(UniqueDate)
            for j=1:length(UniqueDate)
                Urows = find(strcmp(UniqueDate{j},SEDsorted{i,1}(:,InspectionYearCol)));
                SEDsorted{i,1}(Urows,EndIndex2+1)= num2cell(j*ones(length(Urows),1));
                v1=str2num(cell2mat(SEDsorted{i,1}(Urows(1),InspectionYearCol)));
                v2=str2num(cell2mat(SEDsorted{i,1}(Urows(1),ConstructionDateCol)));
                Age=v1-v2;
                SEDsorted{i,1}(Urows,EndIndex2+2)= num2cell(Age*ones(length(Urows),1));
            end
        else
            SEDsorted{i,1}(:,EndIndex2+1)= num2cell(ones(length(SEDsorted{i,1}(:,1)),1));
            if ~isempty(SEDsorted{i,1})
                v1=str2num(cell2mat(SEDsorted{i,1}(1,InspectionYearCol)));
                v2=str2num(cell2mat(SEDsorted{i,1}(1,ConstructionDateCol)));
                Age=v1-v2;
                SEDsorted{i,1}(:,EndIndex2+2)= num2cell(Age*ones(length(SEDsorted{i,1}(:,1)),1));
            end
        end
        %% Unique elements
        if ~isempty(SEDsorted{i,1}(:,[StructuralElemNumSpans,StructuralElemID]))
            [~,ElemID]=uniqueRowsCA(SEDsorted{i,1}(:,[StructuralElemNumSpans,StructuralElemID]));
            jcount=0;
            for j=1:length(ElemID)
                V1=cell2mat(SEDsorted{i,1}(ElemID(j),[StructuralElemNumSpans StructuralElemID]));
                V2=cellstr([char(SEDsorted{i,1}(:,StructuralElemNumSpans)) char(SEDsorted{i,1}(:,StructuralElemID))]);
                V2=strrep(V2,' ','');
                try
                    V2=cellfun(@str2num,V2);
                    EID=find(str2double(V1)==V2);
                catch
                    V1=strrep(V1,' ','');
                    EID=find(strcmp(V1,V2));
                end
%                 V2=cellfun(@str2num,V2);
%                 EID=find(str2double(V1)==V2);
                if ~isempty(EID)
                    jcount=jcount+1;
                    ESEDsorted{i,jcount}=SEDsorted{i,1}(EID,:);
                end
            end
        end
    end
    
    %%
    NTime=2;P2TData={};DeltaData={};
    TransID1=1;TransID2=2; % time steps where the change is calculated
    CID=sum(~cellfun(@isempty,ESEDsorted),2);
    for i=1:length(CID)
        for j=1:CID(i)
            if size(ESEDsorted{i,j},1)>=NTime
                [~,IndSort]=sort(cell2mat(ESEDsorted{i,j}(:,EndIndex2+2)));
                if NTime==1
                    P2TData=[P2TData;ESEDsorted{i,j}(IndSort(1:NTime),:)];
                elseif NTime==2
                    if size(ESEDsorted{i,j},1)>TransID1
                        A1=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID1),EndIndex-4)));
                        A2=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID2),EndIndex-4)));
                        B1=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID1),EndIndex-3)));
                        B2=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID2),EndIndex-3)));
                        C1=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID1),EndIndex-2)));
                        C2=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID2),EndIndex-2)));
                        D1=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID1),EndIndex-1)));
                        D2=str2num(cell2mat(ESEDsorted{i,j}(IndSort(TransID2),EndIndex-1)));
                        DeltaA=A2-A1;
                        DeltaB=B2-B1;
                        DeltaC=C2-C1;
                        DeltaD=D2-D1;
                        M1=0.25*D1+0.5*C1+0.75*B1+A1;
                        M2=0.25*D2+0.5*C2+0.75*B2+A2;
                        DeltaAge=cell2mat(ESEDsorted{i,j}(IndSort(TransID2),EndIndex2+2))-cell2mat(ESEDsorted{i,j}(IndSort(TransID1),EndIndex2+2));
                        ESEDsorted{i,j}(IndSort(TransID1),EndIndex2+3)=ESEDsorted{i,j}(IndSort(TransID1),EndIndex2+2);
                        DeltaData=ESEDsorted{i,j}(IndSort(TransID1),:);
                        DeltaData(1,EndIndex-4)={DeltaA};
                        DeltaData(1,EndIndex-3)={DeltaB};
                        DeltaData(1,EndIndex-2)={DeltaC};
                        DeltaData(1,EndIndex-1)={DeltaD};
                        DeltaData(1,EndIndex2+2)={DeltaAge};
                        DeltaM=M2-M1;
                        DeltaData(1,EndIndex2+4)={DeltaM};
                        DeltaData(1,EndIndex2+5)={M1};
                        DeltaData(1,EndIndex2+6)=ESEDsorted{i,j}(IndSort(TransID2),8); %inspector 2
                        DeltaData(1,EndIndex2+7)={j};
                        P2TData=[P2TData;DeltaData];
                    end
                end
            end
        end
    end
    ESEDsorted(all(cellfun('isempty',ESEDsorted),2),:) = [];
    % StrToRemove={'01869N'};
    % for j=1:length(StrToRemove)
    %      CID=sum(~cellfun(@isempty,ESEDsorted),2);
    %     for i=1:length(CID)
    %         if ~cellfun('isempty',ESEDsorted(i,1))
    %             if strcmp(StrToRemove{j,1},ESEDsorted{i,1}(1,1))
    %                 ESEDsorted(i,:)=[];
    %                 break;
    %             end
    %         end
    %     end
    % end
    if STAtt1==EndIndex-1; STAtt1=[]; end
    if STAtt2==EndIndex-1; STAtt2=[]; end
    if STAtt3==EndIndex-1; STAtt3=[]; end
    if STAtt4==EndIndex-1; STAtt4=[]; end
    if STAtt5==EndIndex-1; STAtt5=[]; end
    if STAtt6==EndIndex-1; STAtt6=[]; end
    if STAtt7==EndIndex-1; STAtt7=[]; end
    if STAtt8==EndIndex-1; STAtt8=[]; end
    if STAtt9==EndIndex-1; STAtt9=[]; end
    
    CID=sum(~cellfun(@isempty,ESEDsorted),2);
%     load('StructuresIDs.mat');
    for i=1:length(CID)
        for j=1:CID(i)
            if ~cellfun('isempty',ESEDsorted(i,1))
                A1=str2double(ESEDsorted{i,j}(:,EndIndex-4));
                B1=str2double(ESEDsorted{i,j}(:,EndIndex-3));
                C1=str2double(ESEDsorted{i,j}(:,EndIndex-2));
                D1=str2double(ESEDsorted{i,j}(:,EndIndex-1));
                M=A1+0.75*B1+0.5*C1+0.25*D1;
                Insp=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,EndIndex)));
                MaterialInd=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,MaterialData)));
                SAtt1=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt1)));
                SAtt2=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt2)));
                SAtt3=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt3)));
                SAtt4=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt4)));
                SAtt5=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt5)));
                SAtt6=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt6)));
                SAtt7=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt7)));
                SAtt8=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt8)));
                SAtt9=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt9)));
                %             if isnan(sum(SAtt1)); SAtt1=[]; end
                %             if isnan(sum(SAtt2)); SAtt2=[]; end
                %             if isnan(sum(SAtt3)); SAtt3=[]; end
                %             if isnan(sum(SAtt4)); SAtt4=[]; end
                %             if isnan(sum(SAtt5)); SAtt5=[]; end
                %             if isnan(sum(SAtt6)); SAtt6=[]; end
                %             if isnan(sum(SAtt7)); SAtt7=[]; end
                %             if isnan(sum(SAtt8)); SAtt8=[]; end
                %             if isnan(sum(SAtt9)); SAtt9=[]; end
                Yr=str2num(cell2mat(ESEDsorted{i,j}(:,2)));
                AGE=Yr-str2num(cell2mat(ESEDsorted{i,j}(:,ConstructionDateCol)));
                [~,iv]=sort(Yr);
                if InterventionCond==3
                    %                 IntCatYr=ones(length(iv),1);
                    %                 InterventionIndex=find(strcmp(string(ESEDsorted{i,j}(1,1)),(table2cell(UISEData(:,1)))));
                    IntYr=zeros(length(iv),1);
                    IntCat=zeros(length(iv),1);
                    InterventionIndex=find(strcmp(string(ESEDsorted{i,j}(1,1)),(table2cell(UISEData(:,1)))));
                    if isempty(InterventionIndex)
                        InterventionIndex=find(str2double(ESEDsorted{i,j}(1,1))==str2double(table2cell(UISEData(:,1))));
                    end
                    if length(InterventionIndex)<j && ~isempty(InterventionIndex)
                        j=1; % in case there are multiple itervetions on one structural element
                    end
                    InterventionYear=str2double(cell2mat(table2cell(UISEData(InterventionIndex(j),2))));
                    InterventionCat=str2double(cell2mat(table2cell(UISEData(InterventionIndex(j),3))));
                    for ii=1:length(InterventionYear)
                        IntIndexTime=find(InterventionYear(ii)>Yr(iv(1)));
                        if ~isempty(IntIndexTime)
                            %                         InterventionTime(ii)=IntIndexTime(ii);
                            %                         IntCatYr(1:InterventionTime(ii))=InterventionYear(ii);
                            %                         IntCatYr(InterventionTime(ii))=InterventionCat(ii);
                            NumObsAfterIntervention(i,j)=length(IntIndexTime);
                            IntYr(1:length(InterventionYear))=InterventionYear;
                            IntCat(1:length(InterventionCat))=InterventionCat;
                        end
                    end
                else
                    InterventionTime=[];
                    IntCatYr=[];
                end
                StrucIndexTrans=find(strcmp(string(ESEDsorted{i,j}(1,1)),string(NewID(:,2))))*ones(length(iv),1);
                if ~isempty(StrucIndexTrans) && InterventionCond>2
                    if sum(strcmp(NoIDIntRows,string(ESEDsorted{i,j}(1,1))))>0
                        if length(InterventionIndex)>1
                            if max(diff(M(iv)))>5
                                IntCat(1)=max(str2double(table2array(UISEData(InterventionIndex,3))));
                            end
                        end
                    end
                    FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) SAtt1 SAtt2 SAtt3 SAtt4 SAtt5 SAtt6 SAtt7 SAtt8 SAtt9 IntYr IntCat];
                else
                    FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) SAtt1 SAtt2 SAtt3 SAtt4 SAtt5 SAtt6 SAtt7 SAtt8 SAtt9];
                end
                %             StrucIndexTrans=find(strcmp(string(ESEDsorted{i,j}(1,1)),string(NewID(:,2))))*ones(length(iv),1);
                %             if ~isempty(StrucIndexTrans)
                %                 FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) SAtt1 SAtt2 SAtt3 SAtt4 SAtt5 SAtt6 SAtt7 SAtt8 SAtt9 IntCatYr];
                %             end
                SSPDsorted{i,j}=FM;
            end
        end
        NetworkB=find(strcmp(string(ESEDsorted{i,j}(1,1)),NetworkData(:,1)));
        NetworkData{NetworkB,ElementType+2}=SSPDsorted(i,:);
    end
    ESEDsorted={};
    SSPDsorted={};
    SEDsorted={};
end

AttributesLabels={app.SelectMaterial.Value; 'Age'; app.SelectSAtt1.Value; app.SelectSAtt2.Value;
    app.SelectSAtt3.Value; app.SelectSAtt4.Value; app.SelectSAtt5.Value; app.SelectSAtt6.Value;
    app.SelectSAtt7.Value; app.SelectSAtt8.Value; app.SelectSAtt9.Value};


[filename, pathname] = uiputfile({'*.mat';'*.*'}, 'Save as','NetworkData');
save(fullfile(pathname,filename),'NetworkData');

filename='Elements';
% [filename, pathname] = uiputfile({'*.mat';'*.*'}, 'Save as',sprintf('ElementsandStructuresData_%s',erase(ColsVal,"/")));
save(fullfile(pathname,filename),'PrimaryElements');

filename='NetworkQuantity';
save(fullfile(pathname,filename),'NetworkData_ElementDetails');

filename='StructrualAttributesLabels';
% [filename, pathname] = uiputfile({'*.mat';'*.*'}, 'Save as',sprintf('StructrualAttributesLabels_%s',erase(ColsVal,"/")));
save(fullfile(pathname,filename),'AttributesLabels');

filename='StructuresID';
% [filename, pathname] = uiputfile({'*.mat';'*.*'}, 'Save as',sprintf('StructrualAttributesLabels_%s',erase(ColsVal,"/")));
save(fullfile(pathname,filename),'NewID');
% if InterventionCond>2
%     filename=sprintf('NIAI_%s',erase(ColsVal,"/"));
%     %     [filename, pathname] = uiputfile({'*.mat';'*.*'}, 'Save as',sprintf('NIAI_%s',erase(ColsVal,"/")));
%     save(fullfile(pathname,filename),'NumObsAfterIntervention');
% end