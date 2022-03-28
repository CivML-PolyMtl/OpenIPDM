ElemAnalyses=0;
Ncurve = app.curve_param;
d.Value = 0; 
YearsDuration=app.TotalYearsDuration;
SelectedYear=app.YearSpinner;
YD=length(YearsDuration);
NumElms=length(Elms);
if NumElms~=0
    EstCat=nan(NumElms,YD);
else
    EstCat=nan(1,YD);
end
for i=1:NumElms
    if ~NetCatAnalyses
        d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
            'Message',sprintf('Processing element %d/%d',i,length(Elms)));
        d.Value = (i/NumElms)/1.2;
        pause(0.001);
    end
    Elm = Elms(i);
    SlCatAnalyses=0;
    [~,IntervType]=RunElementwiseEst(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SelectedYear);
    EstCat(i,1:length(IntCondServiceLife))=IntCondServiceLife;
end
EstCat(EstCat==0)=nan;
ExSL= nanmean(EstCat,1);

if strcmp(app.AnalysesLevelDropDown.Value,'Local') 
    if ~isnan(ExSL)
        cla(app.CatSpeed);
        xticks(app.CatSpeed,'auto');
        yticks(app.CatSpeed,'auto');
        IndServiceLife=(diff(ExSL)>=10^-4);
        XticksVals=1:length(ExSL);
        xlim(app.CatSpeed,[min(XticksVals(IndServiceLife)) max(XticksVals(IndServiceLife))]);
        ylabel(app.CatSpeed,'CDF for Return (Condition)'); 
        plot(app.CatSpeed,XticksVals(IndServiceLife),ExSL(IndServiceLife));
    else
        msgbox('Analyses can not be performed for the selected intervention type.', 'Analyses can not be performed','warn');
    end
end

function [LifeSpanInt,IntervType]=RunElementwiseEst(Elm,app,ElemAnalyses,...
    CatAnalyses,StrucAnalyses,NetAnalyses,SelectedYear)
    SlCatAnalyses=0;
    ButtonState=app.IntServiceLifeButton.Enable;
    app.IntServiceLifeButton.Enable='off';
    da_element();
    app.IntServiceLifeButton.Enable=ButtonState;
end