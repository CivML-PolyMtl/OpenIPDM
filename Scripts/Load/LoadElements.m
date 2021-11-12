FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Elements Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Elements Database file',app.CurrentPath);
end
try
    set(app.InspectionsDatabaseEditField,'value',DataFileName)
    app.SE=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.SE.VariableNames';
    Init_val=[2,3,4,6,7,8,9,10,13,14,15,16,17,19,20];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.SE.SelectedVariableNames = DataColumns(s)';
        app.SE.SelectedFormats(~strcmp(app.SE.SelectedFormats,'%q'))={'%q'};
    end
    app.CurrentPath=FilePath;
catch
    set(app.InspectionsDatabaseEditField,'value','No Selection')
end