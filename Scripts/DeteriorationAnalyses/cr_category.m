ElemAnalyses=0;
Ncurve=app.curve_param;
d.Value = 0; 
YearsDuration=120;
YD=length(YearsDuration);
NumElms=length(Elms);
if NumElms~=0
    CondCr=nan(NumElms,YD);
else
    CondCr=nan(1,YD);
end
First_IndYear = 0;
for i=1:NumElms
    if ~NetCatAnalyses
        d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
            'Message',sprintf('Processing element %d/%d',i,length(Elms)));
        d.Value = (i/NumElms)/1.2;
        pause(0.001);
    end
    Elm = Elms(i);
    [CrCondElem, FirstYear, threshold_value] = RunElementwiseCr(Elm, app, ElemAnalyses, CatAnalyses, StrucAnalyses, NetAnalyses);
    CondCr(i,1:length(CrCondElem)) = CrCondElem;
    First_IndYear = max(First_IndYear,FirstYear);
end
CondCr(CondCr==0)=nan;
if size(CondCr,1) > 1
    ExCr= nanmean(CondCr,1);
else
    ExCr= CondCr;
end
d.close();

ThresholdPDF = diff(ExCr);
IndvYear = find(ExCr >= 0.5);
if isempty(IndvYear)
    true_year = 0;
else
    true_year = 1;
end
Ind_duration = (ThresholdPDF>=10^-4);
Ind_duration = [Ind_duration(1) Ind_duration];
plotting_range = find(Ind_duration);
cla(app.CatSpeed);
if ~isempty(plotting_range)
    xticks(app.CatSpeed,'auto');
    yticks(app.CatSpeed,'auto');
    ylim(app.CatSpeed,[0,1])
    if true_year
        IndThresholdYear = max([IndvYear(1) - plotting_range(1) + 1,1]);
    end
    x_threshold=First_IndYear:First_IndYear+length(ExCr)-1;
    x_threshold=x_threshold(Ind_duration);
    if ~isempty(x_threshold)
        xlim(app.CatSpeed,[x_threshold(1) x_threshold(end)]);
        plot(app.CatSpeed,x_threshold,ExCr(Ind_duration));
        % the index should be shifted becasue diff removes the
        % first year
        if true_year && IndThresholdYear ~= 1
            IndMaxCond = IndThresholdYear ;
            hold(app.CatSpeed)
            plot(app.CatSpeed,[x_threshold(IndMaxCond) x_threshold(IndMaxCond)],[min(ylim(app.CatSpeed)) max(ylim(app.CatSpeed))],'-.');
            text(app.CatSpeed,x_threshold(IndMaxCond)-2, max(ylim(app.CatSpeed))-0.9*max(ylim(app.CatSpeed)), sprintf('$\\tilde{\\mu}_{t|T}= %d \\rightarrow %d$',threshold_value,x_threshold(IndMaxCond)),'Rotation',90,'interpreter','latex','FontSize', 14);
        end
        ylabel(app.CatSpeed,sprintf('Pr[Condition $\\le$ %d]', threshold_value),'interpreter','latex');
        hold(app.CatSpeed)
    else
        msgbox('Analysis is not feasible for this condition threshold.', 'Analyses can not be performed','warn');
    end
else
    msgbox('Analysis is not feasible for this condition threshold.', 'Analyses can not be performed','warn');
end


function [ThresholdCDF, FirstYear, threshold_value]=RunElementwiseCr(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses)
    da_element();
    FirstYear = YearTotal(last_obs_ind);
end