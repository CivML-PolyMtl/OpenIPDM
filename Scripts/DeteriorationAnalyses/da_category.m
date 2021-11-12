ElemAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=app.TotalYearsDuration;
YD=length(YearsDuration);
NumElms=length(Elms);
if NumElms~=0
    ExElm=nan(3,YD,NumElms);
    VarElm=nan(3,3,YD,NumElms);
    YearsElm=nan(NumElms,YD);
    YObs=nan(NumElms,YD);
    Rvals=nan(NumElms,YD);
    QuantityObs=nan(NumElms,YD);
    QuantityModel=nan(NumElms,YD);
    StartYear=max(YearsDuration).*ones(NumElms,1);
else
    ExElm=nan(3,YD,1);
    VarElm=nan(3,3,YD,1);
    YearsElm=nan(1,YD);
    YObs=nan(1,YD);
    Rvals=nan(1,YD);
    QuantityObs=nan(1,YD);
    QuantityModel=nan(1,YD);
    StartYear=max(YearsDuration);
end
for i=1:NumElms
    if CatAnalyses && ~CatReport
        d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
            'Message',sprintf('Processing element %d/%d',i,length(Elms)));
        d.Value = (i/NumElms)/1.2;
        pause(0.001);
    end
    Elm = Elms(i);
    [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Qte]=RunElementwise(Elm,app,...
        ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,CatReport);
    TimeID=find(YearsDuration==YearTotal(1));
    if sum(~isnan(y_Data))>0
        StartYear(i)=YearTotal(1);
    else
        StartYear(i)=max(YearsDuration);
    end
    TimeEnd=length(YearTotal)+TimeID-1;
    ExElm(:,TimeID:TimeEnd,i)=Exsmooth(1:3,:);
    VarElm(:,:,TimeID:TimeEnd,i)=Vsmooth(1:3,1:3,:);
    YearsElm(i,TimeID:TimeEnd)=YearTotal;
    YObs(i,TimeID:TimeEnd)=y_Data;
    Rvals(i,TimeID:TimeEnd)=Re;
    QID=find(~isnan(YObs(i,:)));
    QuantityObs(i,QID)=Qte;
    QuantityModel(i,TimeID:TimeEnd)=Qte;
end
MisStart=min(StartYear);
if NumElms>1
    LambdaModel=QuantityModel./nansum(QuantityModel);
    LambdaObs=QuantityObs./nansum(QuantityObs);
else
    LambdaModel=QuantityModel./(QuantityModel);
    LambdaObs=QuantityObs./(QuantityObs);
end
[Exsmooth,Vsmooth,y_Data,Re]=GaussianMixture_fun(LambdaModel,LambdaObs,...
    ExElm,VarElm,YObs,Rvals);
if CatAnalyses && sum(y_Data~=0)>0
    d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message','Finalizing structural category analyses...');
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
    NotZero=min(find(Exsmooth(1,:)==0))-1;
    Exsmooth=Exsmooth(:,1:NotZero);
    Vsmooth=Vsmooth(:,:,1:NotZero);
    y_Data=y_Data(1:NotZero);
    Re=Re(1:NotZero);
    for j=1:length(Vsmooth)
        Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
    end
    close(d);
    ColorCode=2;
    ylabel(app.CatSpeed,'Speed');
    [xtb,Std,yOr,Rtop,Rlow,x_true]=BackTransformResults(y_Data,Re,Exsmooth,Std,Ncurve,[],100,25);
    InspectorIDLabel_y=round(yOr(~isnan(yOr)));                         % assigned just to have a value 
    PlotTimeSeries(YearTotal,xtb,Std,yOr,Rtop,Rlow,x_true,InspectorIDLabel_y,InterventionVector,app.CatCond,app.CatSpeed,ColorCode);
end
clc
function [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Qte]=RunElementwise(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,CatReport)
    SlCatAnalyses=0; 
    da_element();
end