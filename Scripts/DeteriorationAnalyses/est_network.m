ElemAnalyses=0;
CatAnalyses=0;
StrucAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=app.TotalYearsDuration;
SelectedYear=find(app.YearSpinner.Value==YearsDuration);
YD=length(YearsDuration);
NumStruc=length(Strucs);
EstCat=nan(NumStruc,2);
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
            [ExEst]=RunStructurewiseEst(app,event,ElemAnalyses,...
                CatAnalyses,StrucAnalyses,NetAnalyses,Elms,SelectedYear);
        else
            ExEst=[nan nan];
        end
    else
        ExEst=[nan nan];
    end
    EstCat(i,1:length(ExEst))=round(ExEst,2);
end
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
    'Message','Finalizing network analyses...');
d.Value=0.99;
pause(0.001);
close(d);
app.StructureListBox.Value=CurrentStruc;
app.StructureListBox.ValueChangedFcn(app, event);
app.CategoryListBox.Value=CurrentCat;
NetData=table(Strucs',EstCat(:,1),EstCat(:,2),'VariableNames',{'Bridge ID','Expected Condition','Expected Speed'});
[filename, pathname] = uiputfile({'*.csv';'*.*'}, 'Save as',sprintf('NetworkResults_%s_%d',CurrentCat,app.YearSpinner.Value));
writetable(NetData,fullfile(pathname,filename));
clc
    
function [ExEst]=RunStructurewiseEst(app,event,ElemAnalyses,CatAnalyses,...
    StrucAnalyses,NetAnalyses,Elms,SelectedYear)
    SlCatAnalyses=1;
    NetCatAnalyses=1;
    CatReport=1;
    Elms = app.ElementListBox.Items;
    da_category();
    for j=1:length(Vsmooth)
        Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
    end
    [xtb,Std,yOr,Rtop,Rlow,x_true]=BackTransformResults(y_Data,Re,Exsmooth,...
        Std,Ncurve,[],100,25);[xtb,Std,yOr,Rtop,Rlow,x_true]=BackTransformResults(y_Data,...
        Re,Exsmooth,Std,Ncurve,[],100,25);
    if ~isnan(xtb(2,SelectedYear))
        ExEst=[xtb(1,SelectedYear) xtb(2,SelectedYear)];
    else
        ExEst=[nan nan];
    end
end