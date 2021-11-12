FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Inspectors Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Inspectors Database file',app.CurrentPath);
end
try
    set(app.InspectorsDatabaseEditField,'value',DataFileName)
    app.IP=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.IP.VariableNames';
    Init_val=[2,6];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.IP.SelectedVariableNames = DataColumns(s)';
        app.IP.SelectedFormats(~strcmp(app.IP.SelectedFormats,'%q'))={'%q'};
    end
    app.CurrentPath=FilePath;
catch
    set(app.InspectorsDatabaseEditField,'value','No Selection')
end