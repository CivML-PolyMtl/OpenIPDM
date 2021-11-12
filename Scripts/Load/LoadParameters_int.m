FilterSuffix={'*.mat';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Parameters file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Parameters file',app.CurrentPath);
end
try
    set(app.InterventionParametersEditField,'value',DataFileName)
    app.MPI=load(fullfile(FilePath,DataFileName));
    app.CurrentPath=FilePath;
catch
    set(app.InterventionParametersEditField,'value','No Selection')
end