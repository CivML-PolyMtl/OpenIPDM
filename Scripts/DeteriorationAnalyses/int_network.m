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
                [IntData, Att]=RunStructurewise(app,event,ElemAnalyses,...
                    CatAnalyses,StrucAnalyses,NetAnalyses,Elms);
            for k=1:size(IntData,1)
                if ~isnan(IntData(k,3)) && nansum(IntData(k,:))>0
                    ik = ik + 1;
                    if ik>1
                        Struc_ind_id = find(strcmp(Strucs{i},IntDataStruc(:,1)));
                        year_ind_int = find(IntData(k,4) == cell2mat(IntDataStruc(Struc_ind_id,6)));
                        if isempty(Struc_ind_id) || isempty(year_ind_int)
                            principal_group = 1; secondary_group = 2;
                            [bridge_state_p] = AnalyzeStruc(app, event, IntData(k,4), principal_group);
                            [bridge_state_s] = AnalyzeStruc(app, event, IntData(k,4), secondary_group);
                        else
                            bridge_state_p = [IntDataStruc{Struc_ind_id(end),15}, IntDataStruc{Struc_ind_id(end),16}, IntDataStruc{Struc_ind_id(end),17}, IntDataStruc{Struc_ind_id(end),18},...
                                                IntDataStruc{Struc_ind_id(end),19}, IntDataStruc{Struc_ind_id(end),20}, IntDataStruc{Struc_ind_id(end),21}, IntDataStruc{Struc_ind_id(end),22}];
                            bridge_state_s = [IntDataStruc{Struc_ind_id(end),23}, IntDataStruc{Struc_ind_id(end),24}, IntDataStruc{Struc_ind_id(end),25}, IntDataStruc{Struc_ind_id(end),26}...
                                            IntDataStruc{Struc_ind_id(end),27}, IntDataStruc{Struc_ind_id(end),28}, IntDataStruc{Struc_ind_id(end),29}, IntDataStruc{Struc_ind_id(end),30}];
                        end
                    else
                        principal_group = 1; secondary_group = 2;
                        [bridge_state_p] = AnalyzeStruc(app, event, IntData(k,4), principal_group);
                        [bridge_state_s] = AnalyzeStruc(app, event, IntData(k,4), secondary_group);
                    end
                    IntDataStruc(ik,:)={Strucs{i}, app.CategoryListBox.Items{j},... 
                        IntData(k,1), IntData(k,2), IntData(k,3), IntData(k,4), IntData(k,5),... 
                        IntData(k,6), Att{3}, Att{6}, Att{11}, Att{13}, Att{14}, Att{15},... 
                        bridge_state_p(1), bridge_state_p(2), bridge_state_p(3), bridge_state_p(4),... 
                        bridge_state_p(5), bridge_state_p(6), bridge_state_p(7), bridge_state_p(8),...
                        bridge_state_s(1), bridge_state_s(2), bridge_state_s(3), bridge_state_s(4),...
                        bridge_state_s(5), bridge_state_s(6), bridge_state_s(7), bridge_state_s(8)};
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

if isempty(IntDataStruc)
    IntDataStruc = {nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan,nan};
end
NetIntData=cell2table(IntDataStruc,'VariableNames',{'Bridge ID',...
    'Structural Category','E. Condition','E. Speed','CEC','Int Year','Quantity',...
    'Int Type', 'Municipality', 'Length','Const. Year','DJMA','P. Camion','NumLane',...
    'E.B. Condition (P)','V.B. Condition (P)','E.B. Speed (P)','V.B. Speed (P)',...
    'Diff.B. Condition (P)','Diff.V.B. Condition (P)','Diff.B. Speed (P)','Diff.V.B. Speed (P)',...
    'E.B. Condition (S)','V.B. Condition (S)','E.B. Speed (S)','V.B. Speed (S)','Diff.B. Condition (S)',...
    'Diff.V.B. Condition (S)','Diff.B. Speed (S)','Diff.V.B. Speed (S)' });
[filename, pathname] = uiputfile({'*.csv';'*.*'}, 'Save as',sprintf('NetworkResultsInterventionsDatabase_%s',date()));
try
    writetable(NetIntData,fullfile(pathname,filename));
catch
    message = 'File is not saved';
    uialert(app.MainWindow,message,'Warning','Icon','warning');
end
clc

    
function [IntData, Att]=RunStructurewise(app,event,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,Elms)
    SlCatAnalyses=1;
    NetCatAnalyses=1;
    int_category();
end

function [StrucData]=AnalyzeStruc(app, event, int_year, cat_ind)
    app.FilterCat.Value='Group';
    app.FilterCat.ValueChangedFcn(app,event);
    app.SelectCat.Value=app.SelectCat.Items(cat_ind+1);
    app.SelectCat.ValueChangedFcn(app,event);
    % set initial values
    app.IntServiceLifeButton.Enable=0;
    Strucs=app.StructureListBox.Items;
    Cats=app.CategoryListBox.Items;
    StrucAnalyses=1;NetAnalyses=0;CatReport=0;
    currentSlider = app.ForecastYearsSliderStruc.Value;
    app.ForecastYearsSliderStruc.Value = 15;
    
    AutoCorrectVal=app.AutocorrectElem.Value;
    if strcmp(app.OutliersMissingIntAutocorrectSwitch.Value,'On')
        app.AutocorrectElem.Value='On';
    else
        app.AutocorrectElem.Value='Off';
    end
    addpath(sprintf('%s/Scripts/DeteriorationAnalyses/',pwd))
    run('da_structure.m');
    if sum(Std(1,:))~=0
        int_year_ind = find(int_year-1==YearTotal);
        diffCond = Exsmooth(1,int_year_ind+1) - Exsmooth(1,int_year_ind);
        diffSpeed = Exsmooth(2,int_year_ind+1) - Exsmooth(2,int_year_ind);
        diffCondVar = Vsmooth(1,1,int_year_ind+1) + Vsmooth(1,1,int_year_ind);
        diffSpeedVar = Vsmooth(2,2,int_year_ind+1) + Vsmooth(2,2,int_year_ind);
        stCond = Exsmooth(1,int_year_ind);
        stCondVar = Vsmooth(1,1,int_year_ind);
        stSpeed = Exsmooth(2,int_year_ind);
        stSpeedVar =Vsmooth(2,2,int_year_ind);
    else
        diffCond = NaN;
        diffSpeed = NaN; 
        diffCondVar = NaN; 
        diffSpeedVar = NaN; 
        stCond = NaN; 
        stSpeed = NaN; 
        stCondVar = NaN; 
        stSpeedVar = NaN;
    end
    app.FilterCat.Value='All';
    app.FilterCat.ValueChangedFcn(app,event);
    StrucData = single([stCond stCondVar stSpeed stSpeedVar diffCond diffCondVar diffSpeed diffSpeedVar]);
    app.AutocorrectElem.Value=AutoCorrectVal;
    app.ForecastYearsSliderStruc.Value=currentSlider;
end