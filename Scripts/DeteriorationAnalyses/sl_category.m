ElemAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=120;
YD=length(YearsDuration);
NumElms=length(Elms);
if NumElms~=0
    IntSL=nan(NumElms,YD);
else
    IntSL=nan(1,YD);
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
    [~,IntervType]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses);
    value=find(strcmp(app.IntCodeDropDown.Value,app.IntCodeDropDown.Items))+1;
    if value==IntervType
        SlCatAnalyses=1;
        [IntCondServiceLife,~]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses);
    else
        IntCondServiceLife=nan;
    end
    IntSL(i,1:length(IntCondServiceLife))=IntCondServiceLife;
end
IntSL(IntSL==0)=nan;
ExSL= nanmean(IntSL,1);

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

function [LifeSpanInt,IntervType]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses)
    if SlCatAnalyses
        if app.InterventionCDFButton.Enable
            da_element();
        else
            LifeSpanInt=nan;
        end
    else
        da_element();
        LifeSpanInt=nan;
    end
end