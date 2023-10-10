t=0:Tdb;
if ~isempty(ReTrue) && GraphMatrix(1,1)
    t_=Tdb-1;
    PlotKupdate(t_,Tdb,t,y,xtb,MVv,yearly(1)-1,1,yearly(end),1,...
        OptmInsp,RU,Re,InspBU,Bias_org,InspectorLabel,...
        StructureInd,ElementInd,Xtrue,Rtop,Rlow,RtopTrue,RlowTrue,...
        XTrueParam,MVvTrue,ExsmoothTrue,XbtrueTr,XbbtrueTr); % InpecBiase
    if FigureID>0
        MaxDiff=max(abs(diff(y(~isnan(y)))));
        if MaxDiff>15
            PathFile=sprintf('%s/Tools/PredictRealandSynthetic/AbnormalGraphs',pwd);
        else
            PathFile=sprintf('%s/Tools/PredictRealandSynthetic/Graphs',pwd);
        end
        figure(10)
        set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
            'PaperPositionMode','auto');
        FIG=gcf;
        FIG.PaperPositionMode = 'auto';
        FIG_pos = FIG.PaperPosition;
        FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
        set(gca,'FontSize',20);
        grid on
        saveas(gcf,fullfile(PathFile, sprintf('%d_%dCondition.fig',...
            StructureInd,ElementInd)));
        figure(11)
        set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
            'PaperPositionMode','auto');
        FIG=gcf;
        FIG.PaperPositionMode = 'auto';
        FIG_pos = FIG.PaperPosition;
        FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
        set(gca,'FontSize',20);
        grid on
        saveas(gcf,fullfile(PathFile, sprintf('%d_%dSpeed.fig',...
            StructureInd,ElementInd)));
        figure(12)
        set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
            'PaperPositionMode','auto');
        FIG=gcf;
        FIG.PaperPositionMode = 'auto';
        FIG_pos = FIG.PaperPosition;
        FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
        set(gca,'FontSize',20);
        grid on
        saveas(gcf,fullfile(PathFile, sprintf('%d_%dAcc.fig',...
            StructureInd,ElementInd)));
    end
elseif ~is_synthetic %&& GraphMatrix(1,1)
    if app.BatchMode
    t_=Tdb-1;
    PlotKupdateUI(app,t_,Tdb,t,y,xtb,MVv,yearly(1)-1,1,yearly(end),...
        OptmInsp,RU,Re,InspBU,Bias_org,InspectorLabel,StructureInd,...
        ElementInd,Rtop,Rlow,app.BatchMode,ColorLastObs);
    MaxDiff=max(abs(diff(y(~isnan(y)))));
    if MaxDiff>15
        PathFile=sprintf('%s/Tools/PredictRealandSynthetic/AbnormalGraphs',pwd);
    else
        PathFile=sprintf('%s/Tools/PredictRealandSynthetic/Graphs',pwd);
    end
    figure(10)
    set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
        'PaperPositionMode','auto');
    FIG=gcf;
    FIG.PaperPositionMode = 'auto';
    FIG_pos = FIG.PaperPosition;
    FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
    set(gca,'FontSize',20);
    grid on
    saveas(gcf,fullfile(PathFile, sprintf('%d_%dCondition.fig',...
        StructureInd,ElementInd)));
    figure(11)
    set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
        'PaperPositionMode','auto');
    FIG=gcf;
    FIG.PaperPositionMode = 'auto';
    FIG_pos = FIG.PaperPosition;
    FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
    set(gca,'FontSize',20);
    grid on
    saveas(gcf,fullfile(PathFile, sprintf('%d_%dSpeed.fig',...
        StructureInd,ElementInd)));
    figure(12)
    set(gcf,'Units','inches','Position',[0.5 0.5 10 8],...
        'PaperPositionMode','auto');
    FIG=gcf;
    FIG.PaperPositionMode = 'auto';
    FIG_pos = FIG.PaperPosition;
    FIG.PaperSize = [FIG_pos(3) FIG_pos(4)];
    set(gca,'FontSize',20);
    grid on
    saveas(gcf,fullfile(PathFile, sprintf('%d_%dAcc.fig',...
        StructureInd,ElementInd)));
%elseif FullRun
    end
else
%     if isempty(app)
%         app.BatchMode=0;
%     end
%     for t_=1:Tdb-1
%         PlotKupdateUI(app,t_,Tdb,t,y,xtb,MVv,yearly(1)-1,...
%             1,yearly(end),OptmInsp,RU,Re,InspBU,InpecBiase,...
%             InspectorLabel,StructureInd,ElementInd,Rtop,Rlow,app.BatchMode);
%     end
end
