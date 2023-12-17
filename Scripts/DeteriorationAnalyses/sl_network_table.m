ElemAnalyses=0;
CatAnalyses=1;
StrucAnalyses=0;
Ncurve=app.curve_param;
YearsDuration=120;
YD=length(YearsDuration);
NumStruc=length(Strucs);
CurrentStruc=app.StructureListBox.Value;

d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
        'Message','Preparing for interventions analyses...','Cancelable','on'...
        ,'CancelText','Stop');
d.Value = 0; 
app.AnalysesLevelDropDown.Value = 'Network';
intervention_codes = { 'Intervention Code 1000'
                        'Intervention Code 2000'
                        'Intervention Code 3000'};

all_cat_data = size(app.ModelParam.AllElementsParameters,1);

for iii = 1: length(intervention_codes)
    app.IntCodeDropDown.Value = intervention_codes{iii};
    int_type_report = nan(size(app.ModelParam.AllElementsParameters,1),3);
    if d.CancelRequested
        break
    end
    for ii  = 1: all_cat_data
        CurrentCat = app.ModelParam.AllElementsParameters{ii};
        IntSL=nan(NumStruc,YD);
        if d.CancelRequested
            break
        end
        for i=1:NumStruc
            d.Value = ((ii+(iii-1)*3)/(all_cat_data*3))/1.2;
            d.Message=sprintf('Processing %s for category %d/%d in structure %d/%d',intervention_codes{iii},ii,all_cat_data,i,NumStruc);
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
                else
                    IntCondServiceLife=nan;
                end
            else
                IntCondServiceLife=nan;
            end
            IntSL(i,1:length(IntCondServiceLife))=(IntCondServiceLife);
        end
        
        IntSL(IntSL==0)=nan;
        ExSL= nanmean(IntSL,1);
        app.StructureListBox.Value=CurrentStruc;
        app.StructureListBox.ValueChangedFcn(app, event);

        if nansum(ExSL)>0 && sum(diff(ExSL)>=10^-4)>0
            cla(app.CatSpeed);
            xticks(app.CatSpeed,'auto');
            yticks(app.CatSpeed,'auto');
            IndServiceLife=(diff(ExSL)>=10^-4);
            XticksVals=1:length(ExSL);
            xlim(app.CatSpeed,[min(XticksVals(IndServiceLife)) max(XticksVals(IndServiceLife))]);
            ylabel(app.CatSpeed,'CDF for Return (Condition)');
            plot(app.CatSpeed,XticksVals(IndServiceLife),ExSL(IndServiceLife));

            fig_h = figure('visible','off');
            copyobj(app.CatSpeed, fig_h);
            set(fig_h, 'Resize', 'on');
            figurename = [intervention_codes{iii} '-' erase(CurrentCat,'/')];

            if max(ExSL(IndServiceLife))>0.1 && min(ExSL(IndServiceLife))<0.1
                [~, ind_year_10] = (min(abs(ExSL(IndServiceLife)-0.1)));
            else
                ind_year_10 = nan;
            end
            if ~isnan(ind_year_10)
                year_val10 = XticksVals(IndServiceLife);
                year_val10 = year_val10(ind_year_10);
            else
                year_val10 = nan;
            end

            if max(ExSL(IndServiceLife))>0.5 && min(ExSL(IndServiceLife))<0.5
                [~, ind_year_50] = (min(abs(ExSL(IndServiceLife)-0.5)));
            else
                ind_year_50 = nan;
            end
            if ~isnan(ind_year_50)
                year_val50 = XticksVals(IndServiceLife);
                year_val50 = year_val50(ind_year_50);
            else
                year_val50 = nan;
            end

            if max(ExSL(IndServiceLife))>0.9 && min(ExSL(IndServiceLife))<0.9
                [~, ind_year_90] = (min(abs(ExSL(IndServiceLife)-0.9)));
            else
                ind_year_90 = nan;
            end
            if ~isnan(ind_year_90)
                year_val90 = XticksVals(IndServiceLife);
                year_val90 = year_val90(ind_year_90);
            else
                year_val90 = nan;
            end

            int_type_report(ii,:) = [year_val10, year_val50, year_val90];
            exportgraphics(fig_h,sprintf('Reports/%s/%s.png',FOLDER_NAME,figurename));
        end
    end
    results = table(round(int_type_report(:,1), 2), round(int_type_report(:,2), 2), ...
                    round(int_type_report(:,3), 2), ...
                'VariableNames',{'10th percentile','50th percentile','90th percentile'}); 
    results.Properties.RowNames = app.ModelParam.AllElementsParameters(:,1);
    writetable(results,sprintf('Reports/%s/%s.csv',FOLDER_NAME, intervention_codes{iii}), 'WriteRowNames', true,'Encoding','UTF-8');
end
d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
    'Message','Finalizing network analyses...');
d.Value=0.99;
pause(0.001);
close(d);

clc

    
function [ExSL]=RunStructurewise(app,event,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,Elms)
    SlCatAnalyses=1;
    NetCatAnalyses=1;
    sl_category();
end