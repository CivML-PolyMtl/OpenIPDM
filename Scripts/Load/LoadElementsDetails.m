FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Element Details Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Element Details Database file',app.CurrentPath);
end
try
    set(app.ElementsDetailsEditField,'value',DataFileName)
    app.ED=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.ED.VariableNames';
    Init_val=[1,2,3,5,6,8,10,16,20,26];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.ED.SelectedVariableNames = DataColumns(s)';
        app.ED.SelectedFormats(~strcmp(app.ED.SelectedFormats,'%q'))={'%q'};
    end
    app.CurrentPath=FilePath;
catch
    set(app.ElementsDetailsEditField,'value','No Selection');
end