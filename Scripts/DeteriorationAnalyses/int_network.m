ElemAnalyses=0;
CatAnalyses=1;
StrucAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=120;
YD=length(YearsDuration);
NumStruc=length(Strucs);

IntDataStruc={};%nan(NumStruc,YD);

CurrentStruc=app.StructureListBox.Value;
CurrentCat=app.CategoryListBox.Value;
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message',sprintf('Processing structure %d/%d',i,NumStruc),'Cancelable','on'...
        ,'CancelText','Stop');
ik = 0;
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
        for j = 1:length(app.CategoryListBox.Items)
                app.CategoryListBox.Value=app.CategoryListBox.Items{j};
                app.CategoryListBox.ValueChangedFcn(app, event);
                Elms = app.ElementListBox.Items;
                [IntData]=RunStructurewise(app,event,ElemAnalyses,...
                    CatAnalyses,StrucAnalyses,NetAnalyses,Elms);
            for k=1:size(IntData,1)
                if ~isnan(IntData(k,3)) && nansum(IntData(k,:))>0
                    ik = ik + 1;
                    IntDataStruc(ik,:)={Strucs{i} app.CategoryListBox.Value IntData(k,1) IntData(k,2) IntData(k,3) IntData(k,4) IntData(k,5) IntData(k,6)};
                end
            end
        end
    end
end
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
    'Message','Finalizing network analyses...');
d.Value=0.99;
pause(0.001);
close(d);
app.StructureListBox.Value=CurrentStruc;
app.StructureListBox.ValueChangedFcn(app, event);
app.CategoryListBox.Value=CurrentCat;

NetIntData=cell2table(IntDataStruc,'VariableNames',{'Bridge ID','Structural Category','Elem. Condition','Elem. Speed','CEC','Int Year','Quantity','Int Type'});
[filename, pathname] = uiputfile({'*.csv';'*.*'}, 'Save as',sprintf('NetworkResultsInterventionsDatabase_%s',date()));
try
    writetable(NetIntData,fullfile(pathname,filename));
catch
    message = 'File is not saved';
    uialert(app.MainWindow,message,'Warning','Icon','warning');
end
clc

    
function [IntData]=RunStructurewise(app,event,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,Elms)
    SlCatAnalyses=1;
    NetCatAnalyses=1;
    int_category();
end