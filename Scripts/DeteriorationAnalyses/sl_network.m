ElemAnalyses=0;
CatAnalyses=1;
StrucAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=120;
YD=length(YearsDuration);
NumStruc=length(Strucs);
IntSL=nan(NumStruc,YD);
CurrentStruc=app.StructureListBox.Value;
CurrentCat=app.CategoryListBox.Value;
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message',sprintf('Processing structure %d/%d',i,NumStruc),'Cancelable','on'...
        ,'CancelText','Stop');
for i=1:NumStruc
    d.Value = (i/NumStruc)/1.2;
    d.Message=sprintf('Processing structure %d/%d',i,NumStruc);
    pause(0.001);
    if d.CancelRequested
        break
    end
    app.StructureListBox.Value=Strucs(i);
    app.StructureListBox.ValueChangedFcn(app,event);
    if ~isempty(app.CategoryListBox.Items)
        SelectedCat=find(strcmp(app.CategoryListBox.Items,CurrentCat));
        if ~isempty(SelectedCat)
            app.CategoryListBox.Value=CurrentCat;
            app.CategoryListBox.ValueChangedFcn(app, event);
            Elms = app.ElementListBox.Items;
            [IntCondServiceLife]=RunStructurewise(app,event,ElemAnalyses,...
                CatAnalyses,StrucAnalyses,NetAnalyses,Elms);
        end
    else
        IntCondServiceLife=nan;
    end
    IntSL(i,1:length(IntCondServiceLife))=(IntCondServiceLife);
end
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
    'Message','Finalizing network analyses...');
d.Value=0.99;
pause(0.001);
close(d);
IntSL(IntSL==0)=nan;
ExSL= nanmean(IntSL,1);
app.StructureListBox.Value=CurrentStruc;
app.StructureListBox.ValueChangedFcn(app, event);
app.CategoryListBox.Value=CurrentCat;

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
clc

    
function [ExSL]=RunStructurewise(app,event,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,Elms)
    SlCatAnalyses=1;
    NetCatAnalyses=1;
    sl_category();
end