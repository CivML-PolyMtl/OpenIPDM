function InterventionPrep(NumFiles,InterventionYear)
%global SEDatabaseStore
for i=1:NumFiles
    FilterSuffix={'*.csv';'*.*'};
    [IDataFileName,FilePath]=uigetfile(FilterSuffix,sprintf('Select Intervention Database file %d',i));
    % set(handles.FileNamePath,'String',IDataFileName)
    IDatabaseStore=datastore(fullfile(FilePath,IDataFileName));
    DataColumns=IDatabaseStore.VariableNames';
    [s,v] = listdlg('PromptString','Select Strcture ID and Starting Date to import:',...
        'SelectionMode','multiple',...
        'ListString',DataColumns);
    if v
        ListSelectedColumns=['No Filter';DataColumns(s)];
        IDatabaseStore.SelectedVariableNames = DataColumns(s)';
        for j=1:2
            IDatabaseStore.TextscanFormats{1,s(j)}='%q';
        end
    end
    StoredData{1,i}=readall(IDatabaseStore);
    
    if ~isnumeric(StoredData{1,i}{4,2})
        Dates=table(str2double((StoredData{1,i}{:,2})));
        StoredData{1,i}(:,2)=[];
        StoredData{1,i}(:,2)=Dates;
    end
    StoredData{i}.Properties.VariableNames={'NoStruc','InterventionStart'};
    if i>=2
        StoredData{i}=[StoredData{i-1};StoredData{i}];
    end
end
InterventionData=StoredData{i};
StoredData=[];
[~,UIID]=unique(InterventionData(:,1));
UniqueInterventions=table2cell(InterventionData(UIID,:));
for i=1:length(UniqueInterventions)
    RepeatedIndex=find(str2double(UniqueInterventions(:,1))==str2double(UniqueInterventions(i,1)));
    if length(RepeatedIndex)>1
        [~,selectedID]=max(cell2mat(UniqueInterventions(RepeatedIndex,2)));
        UniqueInterventions(RepeatedIndex(selectedID),:)={0};
    end
end
DeletedInd=cellfun(@(x) (x(:,1)==0), UniqueInterventions);
DInd=find(DeletedInd(:,1));
UniqueInterventions(DInd,:)=[];
SelectedYear=find(cell2mat(UniqueInterventions(:,2))>InterventionYear);
UniqueInterventions(SelectedYear,:)=[];
[filename, pathname] = uiputfile({'*.csv';'*.*'}, 'Save as','MergedInterventions');
writecell(UniqueInterventions,fullfile(pathname,filename));
end