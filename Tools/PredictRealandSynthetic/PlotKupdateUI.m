function PlotKupdateUI(app,t_,T,t,y,E_X,s_X,StartDate,Step,EndDate,...
    OptmInsp,RU,Re,InspBU,InpecBiase,InspectorLabel,StructureInd,...
    ElementInd,TopV,LowV,BatchMode,varargin)
args=varargin;
if ~isempty(args)
    ColorLastObs=args{1};
else
    ColorLastObs=0;
end

figure(10)
if t_+1==T & ColorLastObs
    LastObsInd=find(~isnan(y));
    plot(t(LastObsInd(end)),y(LastObsInd(end)),'*','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','r')
    hold('on');
    plot(t(1:LastObsInd(end-1)),y(1:LastObsInd(end-1)),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','b')
else
    plot(t(1:t_+1),y(1:t_+1),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','b')
end
hold('on');
Rv=zeros(length(Re),1);
Rv(find(OptmInsp==1))=RU(1);
Rv(find(OptmInsp==0))=Re(find(OptmInsp==0));
Rv(find(isnan(OptmInsp)))=NaN;
Rv=[NaN;Rv];
InspB=zeros(length(Re),1);
InspB(find(OptmInsp==1))=y(find(OptmInsp==1))-InspBU;
InspB(find(OptmInsp==0))=y(find(OptmInsp==0))-InpecBiase(find(OptmInsp==0));
InspB(find(isnan(OptmInsp)))=NaN;
% InspB=[NaN;InspB];
InspectorLabel=[NaN;InspectorLabel'];
text(t(1:t_+1), InspB(1:t_+1)-Rv(1:t_+1)-2, num2str(InspectorLabel(1:t_+1)),'FontSize',15);

if t_+1==T & ColorLastObs
    errorbar(t(LastObsInd(end)),InspB(LastObsInd(end)),LowV(LastObsInd(end))',TopV(LastObsInd(end))','LineStyle',...
        'none','CapSize',18,'Color','red','Linewidth',1)
    errorbar(t(1:LastObsInd(end-1)),InspB(1:LastObsInd(end-1)),LowV(1:LastObsInd(end-1))',TopV(1:LastObsInd(end-1))','LineStyle',...
        'none','CapSize',18,'Color','blue','Linewidth',1)
else
    errorbar(t(1:t_+1),InspB(1:t_+1),LowV(1:t_+1)',TopV(1:t_+1)',...
        'LineStyle','none','CapSize',18,'Color','blue','Linewidth',1)
end
plot(t(1:t_+1),E_X(1,1:t_+1),'o--','MarkerFaceColor',[1,1,1],...
    'MarkerEdgeColor','r','LineWidth',1)
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(1,1:t_+1)+s_X(1,1:t_+1),...
    fliplr(E_X(1,1:t_+1)-s_X(4,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(1,1:t_+1)+s_X(5,1:t_+1),...
    fliplr(E_X(1,1:t_+1)-s_X(6,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')
xlabel('Time[Year]')
ylabel(sprintf('Conndition of Structure %d - %d',StructureInd,ElementInd))
%         set('xtick',t)
xlim([0,T-1])
ylim([25,100])
xtickangle(45)
if t_+1==T & ColorLastObs
    h=legend('Hidden Inspection','Inspection','$\pm2 \sigma_{Inspector}$',...
        '$\pm2 \sigma_{Inspector}$','Median','$\pm2 \sigma_{Model}$','$\pm \sigma_{Model}$');
else
    h=legend('Inspection','$\pm2 \sigma_{Inspector}$','Median',...
        '$\pm2 \sigma_{Model}$','$\pm \sigma_{Model}$');
end
set(h,'Interpreter','latex')
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
set(h,'Position',[0.25 0.25 0.1 0.1])
hold('off')
grid('on')
box('on')

figure(11)
plot(t(1:t_+1),E_X(2,1:t_+1),'ro--','MarkerFaceColor',[1,1,1])
hold('on')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(2,1:t_+1)+2*s_X(2,1:t_+1),fliplr(E_X(2,1:t_+1)-2*s_X(2,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(2,1:t_+1)+s_X(2,1:t_+1),fliplr(E_X(2,1:t_+1)-s_X(2,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')
%         set('xtick',t)
xlim([0,T-1])
xlabel('Time[Year]')
ylabel(sprintf('Speed: Structure %d - %d',StructureInd,ElementInd))
xtickangle(45)
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
hold('off')
grid('on')
box('on')

figure(12)
plot(t(1:t_+1),E_X(3,1:t_+1),'ro--','MarkerFaceColor',[1,1,1])
hold('on')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(3,1:t_+1)+2*s_X(3,1:t_+1),fliplr(E_X(3,1:t_+1)-2*s_X(3,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(3,1:t_+1)+s_X(3,1:t_+1),fliplr(E_X(3,1:t_+1)-s_X(3,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')
%         set('xtick',t)
xlim([0,T-1])
% ylim([-0.5,0.5])
xlabel('Time[Year]')
ylabel(sprintf('Acc: Structure %d - %d',StructureInd,ElementInd))
xtickangle(45)
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
hold('off')
grid('on')
box('on')

end
