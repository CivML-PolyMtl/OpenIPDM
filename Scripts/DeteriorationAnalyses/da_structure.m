ElemAnalyses=0;
CatAnalyses=0;
app.estimates_storage{1} = [];
app.estimates_storage{2} = [];
app.estimates_storage{3} = [];
Ncurve=app.curve_param;
YearsDuration=app.TotalYearsDuration;
YD=length(YearsDuration);
NumCats=length(Cats);
if NumCats~=0
    ExElm=nan(3,YD,NumCats);
    VarElm=nan(3,3,YD,NumCats);
    YearsElm=nan(NumCats,YD);
    YObs=nan(NumCats,YD);
    Rvals=nan(NumCats,YD);
    NormModel=zeros(NumCats,YD);
    NormObs=ones(NumCats,YD);
    StartYear=max(YearsDuration).*ones(NumCats,1);
else
    ExElm=nan(3,YD,1);
    VarElm=nan(3,3,YD,1);
    YearsElm=nan(1,YD);
    YObs=nan(1,YD);
    Rvals=nan(1,YD);
    NormModel=zeros(1,YD);
    NormObs=ones(1,YD);
    StartYear=max(YearsDuration);
end

if StrucAnalyses
    d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message',sprintf('Processing Category %d/%d',i,NumCats),'Cancelable','on'...
        ,'CancelText','Stop');
end

d.Value = 0; 

for i=1:length(Cats)
    if StrucAnalyses
        d.Value = (i/NumCats)/1.2;
        d.Message=sprintf('Processing Category %d/%d',i,NumCats);
        if d.CancelRequested
            break
        end
        pause(0.001);
    end
    app.CategoryListBox.Value=Cats(i);
    app.CategoryListBox.ValueChangedFcn(app,event);
    [Exsmooth,Vsmooth,MisStart,y_Data,Re]=RunCategorywise(app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses);
    TimeID=find(YearsDuration==MisStart);
    if nansum(y_Data)>0
        StartYear(i)=MisStart;
    end
    ExElm(:,:,i)=Exsmooth(1:3,:);
    VarElm(:,:,:,i)=Vsmooth(1:3,1:3,:);
    YObs(i,:)=y_Data;
    Rvals(i,:)=Re;
    NormModel(i,:)=(~isnan(Exsmooth(1,:))&Exsmooth(1,:)>0);
end
MisStartStruc=min(StartYear);
if NumCats>1
    LambdaObs=NormObs./sum(~isnan(YObs)&(YObs~=0));
    LambdaModel=NormModel./sum(NormModel);
else
    LambdaObs=NormObs./(~isnan(YObs)&(YObs~=0));
    LambdaModel=NormModel./(NormModel);
end

[Exsmooth,Vsmooth,y_Data,Re]=GaussianMixture_fun(LambdaModel,LambdaObs,...
    ExElm,VarElm,YObs,Rvals);
if StrucAnalyses
    d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message','Finalizing structure analyses...');
    d.Value=0.99;
    pause(0.001);
    InterventionVector=zeros(size(Re));
    NanID=find(y_Data==0);
    y_Data(NanID)=nan;
    Re(NanID)=nan;
    NotNan=min(find(~isnan(y_Data)))-1;
    if ~isempty(NotNan)
        Exsmooth=Exsmooth(:,NotNan:end);
        Vsmooth=Vsmooth(:,:,NotNan:end);
        y_Data=y_Data(NotNan:end);
        Re=Re(NotNan:end);
    end
    NotZero=min(find(Exsmooth(1,:)==0))-1;
    if NotZero~=0
        Exsmooth=Exsmooth(:,1:NotZero);
        Vsmooth=Vsmooth(:,:,1:NotZero);
        y_Data=y_Data(1:NotZero);
        Re=Re(1:NotZero);
    end
    for j=1:length(Vsmooth)
        Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
    end
    if sum(Std(1,:))~=0
        YearTotal=MisStartStruc:MisStartStruc+NotZero-1;
        close(d);
        ColorCode=3;
        [xtb,Std,yOr,Rtop,Rlow,y_empty]=BackTransformResults(y_Data,Re,Exsmooth,Std,Ncurve,[],100,25);
        InspectorIDLabel_y=round(yOr(~isnan(yOr)));                      % assigned just to have a value 
        PlotTimeSeries(YearTotal,xtb,Std,yOr,[],Rtop,Rlow,[],InspectorIDLabel_y,InterventionVector,app.StrucCond,app.StrucSpeed,ColorCode);
        app.estimates_storage{1} = YearTotal;
        app.estimates_storage{2} = xtb;
        app.estimates_storage{3} = Std;
    else
        message = 'The selected structure has no inspection data, or it has no model associated with it.';
        uialert(app.MainWindow,message,'Warning','Icon','warning');
        close(d);
        app.estimates_storage{1} = [];
        app.estimates_storage{2} = [];
        app.estimates_storage{3} = [];
        clc
    end
end
clc
function [Exsmooth,Vsmooth,MisStart,y_Data,Re]=RunCategorywise(app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses)
    Elms = app.ElementListBox.Items;
    CatReport=0;
    da_category();
end