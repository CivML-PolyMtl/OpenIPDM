FilterSuffix={'*.csv';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Structures Database file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Structures Database file',app.CurrentPath);
end
try
    set(app.StructuresDatabaseEditField,'value',DataFileName)
    app.ST=datastore(fullfile(FilePath,DataFileName));
    DataColumns=app.ST.VariableNames';
    Init_val=[1,4,7,11,12,13,14,15,16,17,18,19,23,24,25];
    [s,v] = listdlg('PromptString','Select columns to import:',...
        'SelectionMode','multiple','InitialValue',Init_val,...
        'ListString',DataColumns);
    if v
        app.ST.SelectedVariableNames = DataColumns(s)';
        app.ST.SelectedFormats(~strcmp(app.ST.SelectedFormats,'%q'))={'%q'};
    end

    app.CurrentPath=FilePath;
catch
    set(app.StructuresDatabaseEditField,'value','No Selection')
end