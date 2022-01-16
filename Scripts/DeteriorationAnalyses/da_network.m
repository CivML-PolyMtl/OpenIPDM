ElemAnalyses=0;
CatAnalyses=0;
StrucAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=app.TotalYearsDuration;
YD=length(YearsDuration);
NumStruc=length(Strucs);
ExElm=nan(3,YD,NumStruc);
VarElm=nan(3,3,YD,NumStruc);
YearsElm=nan(NumStruc,YD);
YObs=nan(NumStruc,YD);
Rvals=nan(NumStruc,YD);
NormModel=zeros(NumStruc,YD);
NormObs=ones(NumStruc,YD);
 d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message',sprintf('Processing structure %d/%d',i,NumStruc),'Cancelable','on'...
        ,'CancelText','Stop');
for i=1:NumStruc
   
    d.Value = (i/NumStruc)/1.2;
    d.Message=sprintf('Processing structure %d/%d',i,NumStruc);
    if d.CancelRequested
        break
    end
    pause(0.001);
    app.StructureListBox.Value=Strucs(i);
    app.StructureListBox.ValueChangedFcn(app,event);
    if strcmp(app.StructuralElementsGroupNet.Value,'Principal')
        app.FilterCat.Value='Group';
        app.FilterCat.ValueChangedFcn(app,event);
        if sum(strcmp('Principal',app.SelectCat.Items))>1
            app.SelectCat.Value='Principal';
            app.SelectCat.ValueChangedFcn(app,event);
        end
    else
        app.FilterCat.Value='Group';
        app.FilterCat.ValueChangedFcn(app,event);
        if sum(strcmp('Secondaire',app.SelectCat.Items))>1
            app.SelectCat.Value='Secondaire';
            app.SelectCat.ValueChangedFcn(app,event);
        end
    end
    [Exsmooth,Vsmooth,MisStart,y_Data,Re]=RunStructurewise(app,event,ElemAnalyses,...
        CatAnalyses,StrucAnalyses,NetAnalyses);
    TimeID=find(YearsDuration==MisStart);
    MinTimeSeries(i)=min(find(Exsmooth(1,TimeID:end)==0))-1;
    StartYear(i)=MisStart;
    ExElm(:,:,i)=Exsmooth(1:3,:);
    VarElm(:,:,:,i)=Vsmooth(1:3,1:3,:);
    YObs(i,:)=y_Data;
    Rvals(i,:)=Re;
    NormModel(i,:)=(~isnan(Exsmooth(1,:))&Exsmooth(1,:)>0);
end
MisStartStruc=min(StartYear);
MinIndex=min(MinTimeSeries(MinTimeSeries>1));
if NumStruc>1
    LambdaObs=NormObs./sum(~isnan(YObs)&(YObs~=0));
    LambdaModel=NormModel./sum(NormModel);
else
    LambdaObs=NormObs./(~isnan(YObs)&(YObs~=0));
    LambdaModel=NormModel./(NormModel);
end

[Exsmooth,Vsmooth,y_Data,Re]=GaussianMixture_fun(LambdaModel,LambdaObs,...
    ExElm,VarElm,YObs,Rvals);
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
    'Message','Finalizing network analyses...');
d.Value=0.99;
pause(0.001);
InterventionVector=zeros(size(Re));
NanID=find(y_Data==0);
y_Data(NanID)=nan;
Re(NanID)=nan;
NotNan=min(find(~isnan(y_Data)))-1;
Exsmooth=Exsmooth(:,NotNan:end);
Vsmooth=Vsmooth(:,:,NotNan:end);
y_Data=y_Data(NotNan:end);
Re=Re(NotNan:end);
% NotZero=min(find(Exsmooth(1,:)==0))-1;
Exsmooth=Exsmooth(:,1:MinIndex);
Vsmooth=Vsmooth(:,:,1:MinIndex);
y_Data=y_Data(1:MinIndex);
Re=Re(1:MinIndex);
for j=1:length(Vsmooth)
    Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
end
YearTotal=MisStartStruc:MisStartStruc+MinIndex-1;
close(d);
[xtb,Std,yOr,Rtop,Rlow,y_empty]=BackTransformResults(y_Data,Re,Exsmooth,Std,Ncurve,[],100,25);
InspectorIDLabel_y=round(yOr(~isnan(yOr)));
app.NetworkResults.Ex=xtb;
app.NetworkResults.Std=Std;
app.NetworkResults.YearTotal=YearTotal;
app.NetworkResults.yOr=yOr;
app.NetworkResults.Rtop=Rtop;
app.NetworkResults.Rlow=Rlow;
app.NetworkResults.InspectorIDLabel_y=InspectorIDLabel_y;
app.NetworkResults.InterventionVector=InterventionVector;
app.NetworkResults.Status=app.StatusDropDownNetwork.Value;
app.NetworkResults.Group=app.StructuralElementsGroupNet.Value;
app.NetworkResults.FilterStruc=app.FilterStructures.Value;
app.NetworkResults.StrucVal=app.StructValue.Value;
app.NetworkResults.Slider=app.ForecastYearsSliderNet.Value;
ColorCode=4;
PlotTimeSeries(YearTotal,xtb,Std,yOr,[],Rtop,Rlow,[],InspectorIDLabel_y,InterventionVector,app.NetCond,app.NetSpeed,ColorCode);
clc    
function [Exsmooth,Vsmooth,MisStartStruc,y_Data,Re]=RunStructurewise(app,event,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses)
    StrucStatus=unique(app.BridgeAttributes(:,12));
    % Exclude or include a structure according to its status
    if strcmp(StrucStatus,app.StatusDropDownNetwork.Value)
        Cats = app.CategoryListBox.Items;
    else
        Cats=[];
    end
    da_structure();
end