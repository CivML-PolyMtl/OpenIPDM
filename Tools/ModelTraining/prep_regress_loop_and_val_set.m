MultiPass = 20;
RegressLoop=1;
SumLL=zeros(MultiPass,1);
Converge=0;
TotalTimeSteps=size(y,3);
OriginalValues=100;
[~,MaxCondition]=SpaceTransformation(Ncurve,OriginalValues,100,25);
NormObsValues = gather(y(1,:,2));
NormObsValues(find(NormObsValues>MaxCondition(end)))=MaxCondition(end);

y_valid=MdataEngy.ModelValid.YS;
Re_Valid=MdataEngy.ModelValid.ReS;
R_valid=MdataEngy.ModelValid.RS;
RU_valid=MdataEngy.ModelValid.RUS;
InpecBiase_valid=MdataEngy.ModelValid.InpecBiaseS;
InspBU_valid=MdataEngy.ModelValid.InspBUS;
ObsYears_valid=MdataEngy.ModelValid.ObsYearsS;
InspectorsObs_valid=MdataEngy.ModelValid.AllInspectors;
for i=1:length(InspectorsID)
    ReReshape_valid=reshape(Re_Valid,[size(Re_Valid,2),size(Re_Valid,3)]);
    ReReshape_valid(find(InspectorsObs_valid{i}))=EngBiasData(i,end).^2;
    Re_Valid(1,:,:)=ReReshape_valid;
    % update inspector bias
    InpecBiaseReshape_valid=reshape(InpecBiase_valid,[size(InpecBiase_valid,2),size(InpecBiase_valid,3)]);
    InpecBiaseReshape_valid(find(InspectorsObs_valid{i}))=EngBiasData(i,2);
    InpecBiase_valid(1,:,:)=InpecBiaseReshape_valid;
end
CurrentInspectorObs_valid=zeros(size(InspectorsObs_valid{1}),'gpuArray');
init_V_valid=zeros(3,3,length(MdataEngy.ModelValid.RS),'gpuArray');
if ~isempty(CurrentInspectorID)
    RUReshape_valid=reshape(Re_Valid,[size(Re_Valid,2),size(Re_Valid,3)]);
    CurrentInspectorIndex_valid=find(CurrentInspectorID==InspectorsID);
    CurrentInspectorObs=InspectorsObs_valid{CurrentInspectorIndex_valid};
    RUReshape_valid(find(CurrentInspectorObs))=(CurrentInspectorParam(1)).^2;
    RU_valid(1,:,:)=RUReshape_valid;
    % update inspector bias
    InspBUReshape_valid=reshape(InspBU_valid,[size(InspBU_valid,2),size(InspBU_valid,3)]);
    InspBUReshape_valid(find(CurrentInspectorObs_valid))=(CurrentInspectorParam(2));
    InspBU_valid(1,:,:)=InspBUReshape_valid;
end