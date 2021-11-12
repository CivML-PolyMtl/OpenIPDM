FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Interventions Budget Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Interventions Budget Database file',app.CurrentPath);
end
try
    set(app.BudgetDatabaseEditField,'value',DataFileName)
    app.IB=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.IB.VariableNames';
    Init_val=[1,2,3];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.IB.SelectedVariableNames = DataColumns(s)';
        app.IB.SelectedFormats(~strcmp(app.IB.SelectedFormats,'%q'))={'%q'};
    end
    app.CurrentPath=FilePath;
catch
    set(app.BudgetDatabaseEditField,'value','No Selection');
end