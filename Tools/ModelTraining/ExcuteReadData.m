% Database Columns
% 1: Bridge ID % 2: Structure Type % 3: Munipalite % 4: LAT % 5: Long
% 6: Longtotale % 7: LongTablier % 8: LargHorstout % 9: LargCarrossable
% 10: SuperfTablier % 11: AnConst % 12: Statut % 13: DJMA % 14: x_camions
% 15: NBREVoies.
% Deterioration Info:
% 1: NOStruc % 2: DateInsp % 3: NoTrav % 4: Element % 5: NoElem % 6: Position
% 7: Materiau % 8: TypeElement % 9: A % 10: B % 11: C % 12: D % 13: CEC
% 14: Inspection ID (rreplaced by Inspector ID) % 15: Element ID

SV=SE.SelectedVariableNames;
SE.SelectedFormats(~strcmp(SE.SelectedFormats,'%q'))={'%q'};
SEData=readall(SE);
FilterCond=SelectedElementCat;%~= 1

ColsVal=SelectedElementCat;
ColsInd=4; % Element
Rows=find(strcmp(ColsVal,table2array(SEData(:,ColsInd))));
SEData=SEData(Rows,SV);
INData=readall(IN);
InterventionsRows=find(strcmp(ColsVal,table2array(INData(:,ColsInd+1))));
INData=INData(InterventionsRows,[1 2 4 6]);



%% check interventions
if InterventionCond~=1
    AllRows=[];
    InterAllRows=[];
    NoIDIntRows=[];
    UISEData=unique(INData);
    ISID=str2double(table2cell(SEData(:,end)));
    ISID_ind=str2double(table2cell(SEData(:,1)));
    for i=1:size(table2cell(UISEData),1)
        SRows = find(strcmp(table2cell(UISEData(i,end)),table2cell(SEData(:,end))));
        if isempty(SRows)
            SRows=find(ISID==str2double(cell2mat(table2cell(UISEData(i,end)))));
        end
        if isempty(SRows)
            SRows = find(strcmp(table2cell(UISEData(i,1)),table2cell(SEData(:,1))));
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
if ~isempty(SEData)
    MaterialData=7;
    if MaterialData~=0
        MaterialCategories=unique(SEData(:,MaterialData));
        MaterialCategories=table2array(MaterialCategories);
        if ~isnumeric(MaterialCategories(1,1))
            for i=1:length(MaterialCategories)
                SEData(find(strcmp(MaterialCategories{i},table2array(SEData(:,MaterialData)))),MaterialData)={sprintf('%d',i)};
            end
        end
    end
    TypeElementData=8;
    %% Specify Type of Element Index
    if TypeElementData~=0
        TypeElementCategories=unique(SEData(:,TypeElementData));
        TypeElementCategories=table2array(TypeElementCategories);
        if ~isnumeric(TypeElementCategories(1))
            for i=1:length(TypeElementCategories)
                SEData(find(strcmp(TypeElementCategories{i},table2array(SEData(:,TypeElementData)))),TypeElementData)={sprintf('%d',i)};
            end
        end
    end
end
if InterventionCond==2
    NoInterRows=setdiff(1:size(SEData,1),AllRows)';
    SEData=SEData(NoInterRows,:);
elseif InterventionCond==3
    SEData=SEData(AllRows,:);
    UISEData=UISEData(InterAllRows,:);
end
SEData=table2array(SEData);
if ~isempty(SEData)
    %% Bridge Database
    ST.SelectedFormats(strcmp(ST.SelectedFormats,'%f'))={'%q'};
    StructuresData=readall(ST);
    NewID(:,1)=1:height(StructuresData(:,1));
    NewID=num2cell(NewID);
    NewID(:,2)=table2cell(StructuresData(:,1));
    
    if ~strcmp(StatusValue,'All')
        StatusRows=find(strcmp(StatusValue,table2array(StructuresData(:,12))));
        %     StrucVars=app.SelectStructureActiveCol.Items;
        %     StrucVars(1)=[];
        StructuresData=StructuresData(StatusRows,ST.SelectedVariableNames);
    end
    StructuresData=table2array(StructuresData);
    
    %% Inspectors
    InspectorsDatabase=readall(IP);
    InspectorsDatabase=table2array(InspectorsDatabase);
    EmptyEng=[];
    for i=1:size(SEData,1)
        EngId=find(strcmp(SEData(i,end-1),InspectorsDatabase(:,1)));
        if ~isempty(EngId)
            SEData(i,end-1)=num2cell(InspectorsDatabase(EngId,2));
        else
            EmptyEng=[EmptyEng;i];
        end
    end
    %% add construction date to elements databse
    SEData(EmptyEng,:)=[];
    EndIndex=size(SEData,2);
    for i=1:length(StructuresData)
        Ind=find(strcmp(StructuresData(i,1),SEData(:,1)));
        for j=1:size(StructuresData,2)-3
            if ~isempty(Ind)
                SEData(Ind,EndIndex+j)=StructuresData(i,j+3); % Const Date
            end
        end
    end
    
    
    
    InspectionYearCol=2;
    StructuralElemID=5;
    StructuralElemNumSpans=3;
    
    
    
    STAtt1=1+EndIndex;% Lat
    STAtt2=2+EndIndex;% Long
    STAtt3=3+EndIndex;% Longtotale
    STAtt4=4+EndIndex;% LongTablier
    STAtt5=5+EndIndex;% LargHorstout
    STAtt6=6+EndIndex;% LargCarrossable
    STAtt7=7+EndIndex;% SuperfTablier
    ConstructionDateCol=8+EndIndex; % AnConst
    STAtt8=10+EndIndex;% DJMA
    STAtt9=11+EndIndex;% x_camions
    STAtt10=12+EndIndex;% NBREVOIES

    
    MetaData.AttributesLabels={'Material';'Age'; 'Lat'; 'Long';'Longtotale'; 'LongTablier';
        'LargHorstout'; 'LargCarrossable';'SuperfTablier'; 'DJMA'; 'x_camions'; 'NBREVOIES'};
    MetaData.Material=MaterialCategories;
    MetaData.TypeElement=TypeElementCategories;
    
    if InterventionCond==3
        save(sprintf('%s/ExtractedData/MetaData_Intervention_%s.mat',OriginPWD,erase(ColsVal,"/")),'MetaData');
    else
        save(sprintf('%s/ExtractedData/MetaData_%s.mat',OriginPWD,erase(ColsVal,"/")),'MetaData');
    end
    
    %% Convert inspecction Date data
    for i=1:size(SEData,1)
        DateVec=datestr(datenum(SEData(i,InspectionYearCol)));
        SEData(i,InspectionYearCol)={DateVec(end-3:end)};
    end
    
    %% Sort Elements to Structures
    EndIndex2=size(SEData,2);
    for i=1:length(StructuresData(:,1))
        SRows = find(strcmp(StructuresData(i,1),SEData(:,1)));
        SEDsorted{i,1}=SEData(SRows,:);
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
        if ~isempty(SEDsorted{i,1}(:,15))
            [ElementsNum,ElemID]=unique(SEDsorted{i,1}(:,15));
            jcount=0;
            for j=1:length(ElemID)
                V1=cell2mat(SEDsorted{i,1}(ElemID(j),15));
                V2=cellstr(SEDsorted{i,1}(:,15));
                EID=find(strcmp(V1,V2));
                if ~isempty(EID)
                    jcount=jcount+1;
                    ESEDsorted{i,jcount}=SEDsorted{i,1}(EID,:);
                end
            end
            %% added on 4/6/2023
        elseif ~isempty(SEDsorted{i,1}(:,[StructuralElemNumSpans,StructuralElemID]))
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
    ESEDsorted(all(cellfun('isempty',ESEDsorted),2),:) = [];
    
    if STAtt1==EndIndex-1; STAtt1=[]; end
    if STAtt2==EndIndex-1; STAtt2=[]; end
    if STAtt3==EndIndex-1; STAtt3=[]; end
    if STAtt4==EndIndex-1; STAtt4=[]; end
    if STAtt5==EndIndex-1; STAtt5=[]; end
    if STAtt6==EndIndex-1; STAtt6=[]; end
    if STAtt7==EndIndex-1; STAtt7=[]; end
    if STAtt8==EndIndex-1; STAtt8=[]; end
    if STAtt9==EndIndex-1; STAtt9=[]; end
    if STAtt10==EndIndex-1; STAtt10=[]; end
    
    CID=sum(~cellfun(@isempty,ESEDsorted),2);
    
    %load('StructuresIDs.mat');
    for i=1:length(CID)
        for j=1:CID(i)
            if ~cellfun('isempty',ESEDsorted(i,1))
                A1=str2double(ESEDsorted{i,j}(:,EndIndex-6));
                B1=str2double(ESEDsorted{i,j}(:,EndIndex-5));
                C1=str2double(ESEDsorted{i,j}(:,EndIndex-4));
                D1=str2double(ESEDsorted{i,j}(:,EndIndex-3));
                M=A1+0.75*B1+0.5*C1+0.25*D1;
                Insp=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,EndIndex-1)));
                if ~isempty(categorical_data_ind)
                    MaterialInd = str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,categorical_data_ind)));
                end
                SAtt1=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt1)));
                SAtt2=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt2)));
                SAtt3=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt3)));
                if ~sum(strcmp(SelectedElementCat,{'Chasse-roue / trottoir',...
                        'Dessous de la dalle/voûte','Mur','Mur de tête',...
                        'Mur en aile','Murs/naiss.voûte/coins infér.',...
                        'Plafond suspendu - Tuiles','Revêtement de mur',...
                        'Voûte / Dalle'}))
                    SAtt4=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt4)));
                else
                    SAtt4 = [];
                end
                if ~sum(strcmp(SelectedElementCat,{'Dessous de la dalle/voûte','Mur','Mur de tête','Mur en aile','Murs/naiss.voûte/coins infér.'}))
                    SAtt5=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt5)));
                else
                    SAtt5 = [];
                end
                if ~sum(strcmp(SelectedElementCat,{'Mur'}))
                    SAtt6=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt6)));
                else
                    SAtt6 = [];
                end
                SAtt7=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt7)));
                if ~sum(strcmp(SelectedElementCat,{'Mur'}))
                    SAtt8=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt8)));
                else
                    SAtt8 = [];
                end
                if ~sum(strcmp(SelectedElementCat,{'Mur'}))
                    SAtt9=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt9)));
                else
                    SAtt9 = [];
                end
                if ~sum(strcmp(SelectedElementCat,{'Tirants','Toiture'}))
                    SAtt10=str2double(cellfun(@cellstr,ESEDsorted{i,j}(:,STAtt10)));
                else
                    SAtt10 = [];
                end

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
                    if length(InterventionIndex)<j
                        j=1; % in case there are multiple itervetions on one structural element
                    end
                    InterventionYear=str2double(cell2mat(table2cell(UISEData(InterventionIndex(j),2))));
                    InterventionCat=str2double(cell2mat(table2cell(UISEData(InterventionIndex(j),3))));
                    for ii=1:length(InterventionYear)
                        IntIndexTime=find(InterventionYear(ii)<=Yr);
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
                All_SAtt = [SAtt1 SAtt2 SAtt3 SAtt4 SAtt5 SAtt7 SAtt8 SAtt9 SAtt10]; %SAtt6
                
                if ~isempty(StrucIndexTrans) && InterventionCond>2
                    if sum(strcmp(NoIDIntRows,string(ESEDsorted{i,j}(1,1))))>0
                        if length(InterventionIndex)>1
                            if max(diff(M(iv)))>5
                                IntCat(1)=max(str2double(table2array(UISEData(InterventionIndex,3))));
                            end
                        end
                    end
                    FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) All_SAtt IntYr IntCat];
                else
                    FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) All_SAtt];
                end
                %             StrucIndexTrans=find(strcmp(string(ESEDsorted{i,j}(1,1)),string(NewID(:,2))))*ones(length(iv),1);
                %             if ~isempty(StrucIndexTrans)
                %                 FM=[M(iv) Insp(iv) Yr(iv) StrucIndexTrans MaterialInd AGE(iv) SAtt1 SAtt2 SAtt3 SAtt4 SAtt5 SAtt6 SAtt7 SAtt8 SAtt9 IntCatYr];
                %             end
                if IncludeStructuralAtt
                    if sum(any(isnan(All_SAtt)))==0 && sum(any(AGE<0))==0
                        SSPDsorted{i,j}=FM;
                    end
                else
                    SSPDsorted{i,j}=FM;
                end
            end
        end
        if exist('SSPDsorted','var') == 1
            if i<=size(SSPDsorted,1)
                ind_full = find(~cellfun(@isempty,SSPDsorted(i,:)));
            end
            if min(ind_full) ~= 1
                 place_holder = SSPDsorted(i,ind_full);
                 SSPDsorted(i,ind_full) = {[]};
                 SSPDsorted(i,1:length(ind_full)) = place_holder;
            end
        end
    end
    if exist('SSPDsorted','var') == 0
        SSPDsorted = {[]};
    end
    CID=sum(~cellfun(@isempty,SSPDsorted),2);
    empty_ind = find(sum(~cellfun(@isempty,SSPDsorted),2)==0);
    if ~isempty(empty_ind)
        CID(empty_ind)=[];
        SSPDsorted(empty_ind,:)=[];
    end
    Inspectors=[];
    for i=1:length(CID)
        for j=1:CID(i)
            if ~isempty(SSPDsorted{i,j})
                UniqueInspectors = unique(SSPDsorted{i,j}(:,2));
                Inspectors = union(UniqueInspectors,Inspectors); 
            end
        end
    end
    TrainingWindow=50;
    [Inspectors]=FIlteredInspectors(SSPDsorted,Inspectors,TrainingWindow);
    
    if InterventionCond==3
        FullPathEx=sprintf('%s/ExtractedData/InspectionData_Intervention_%s.mat',OriginPWD,erase(ColsVal,"/"));
        FullPath_Inspectors=sprintf('%s/ExtractedData/Inspectors_int_%s.mat',OriginPWD,erase(ColsVal,"/"));
        save(sprintf('%s/ExtractedData/NumAfterIntervention_%s.mat',OriginPWD,erase(ColsVal,"/")),'NumObsAfterIntervention');
        save(FullPathEx,'SSPDsorted');
        save(sprintf('%s/ExtractedData/MetaData_Intervention_%s.mat',OriginPWD,erase(ColsVal,"/")),'MetaData');
        save(FullPath_Inspectors,'Inspectors');
        save(sprintf('%s/ExtractedData/StructuresID_Intervention_%s.mat',OriginPWD, erase(ColsVal,"/")), 'NewID');
    else
        FullPathEx=sprintf('%s/ExtractedData/InspectionData_%s.mat',OriginPWD,erase(ColsVal,"/"));
        FullPath_Inspectors=sprintf('%s/ExtractedData/Inspectors_%s.mat',OriginPWD,erase(ColsVal,"/"));
        save(FullPathEx,'SSPDsorted');
        save(sprintf('%s/ExtractedData/MetaData_%s.mat',OriginPWD,erase(ColsVal,"/")),'MetaData');
        save(FullPath_Inspectors,'Inspectors');
        save(sprintf('%s/ExtractedData/StructuresID_%s.mat',OriginPWD, erase(ColsVal,"/")), 'NewID');
    end
else
    FullPathEx='';
    FullPath_Inspectors='';
end
