FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Interventions Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Interventions Database file',app.CurrentPath);
end
try
    set(app.InterventionsDatabaseEditField,'value',DataFileName)
    app.IN=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.IN.VariableNames';
    Init_val=[2,4,5,8,13,16];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.IN.SelectedVariableNames = DataColumns(s)';
        app.IN.SelectedFormats(~strcmp(app.IN.SelectedFormats,'%q'))={'%q'};
    end
    app.CurrentPath=FilePath;
catch
    set(app.InterventionsDatabaseEditField,'value','No Selection')
end