FilterSuffix={'*.mat';'*.*'};
if isempty(app.CurrentPath)
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Parameters file');
else
    [DataFileName,FilePath]=uigetfile(FilterSuffix,'Select Parameters file',app.CurrentPath);
end
try
    set(app.ModelParametersEditField,'value',DataFileName)
    app.MP=load(fullfile(FilePath,DataFileName));
    app.CurrentPath=FilePath;
catch
    set(app.ModelParametersEditField,'value','No Selection')
end
