function PlotKupdate(t_,T,t,y,E_X,s_X,StartDate,Step,EndDate,Fi,OptmInsp,...
    RU,Re,InspBU,InpecBiase,InspectorLabel,StructureInd,ElementInd,XTrue,...
    TopV,LowV,RtopTrue,RlowTrue,ReOrTrue,MVvTrue,ExsmoothTrue,XbtrueTr,...
    XbbtrueTr)

figure(10)
set(gcf,'units','normalized','outerposition',[0 0 1 1])
set(gcf,'Units','inches','Position',[0.5 0.5 6.5 4.5],'PaperPositionMode','auto');
plot(t(1:t_+1),y(1:t_+1),'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','b')
hold on

Rv=zeros(length(Re),1);
Rv(find(OptmInsp==1))=RU(find(OptmInsp==1));
Rv(find(OptmInsp==0))=Re(find(OptmInsp==0));
Rv(find(isnan(OptmInsp)))=NaN;
Rv=[NaN;Rv];

InspB=zeros(length(Re),1);
InspB(find(OptmInsp==1))=y(find(OptmInsp==1))-InspBU;
InspB(find(OptmInsp==0))=y(find(OptmInsp==0))-InpecBiase(find(OptmInsp==0));
InspB(find(isnan(OptmInsp)))=NaN;

plot(t(1:t_+1),InspB(1:t_+1),'b*')

InspectorLabel=[NaN;InspectorLabel'];
text(t(1:t_+1), InspB(1:t_+1)-Rv(1:t_+1)-2, num2str(InspectorLabel(1:t_+1)));
errorbar(t(1:t_+1),InspB(1:t_+1),LowV(1:t_+1)',TopV(1:t_+1)','LineStyle','none','CapSize',20,'Color','blue','Linewidth',1)
errorbar(t(1:t_+1),InspB(1:t_+1),RlowTrue(1:t_+1)',RtopTrue(1:t_+1)','LineStyle','none','CapSize',10,'Color','black','Linewidth',1)

plot(t(1:t_+1),E_X(1,1:t_+1),'o--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','r','LineWidth',1)
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(1,1:t_+1)+s_X(1,1:t_+1),fliplr(E_X(1,1:t_+1)-s_X(4,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(1,1:t_+1)+s_X(5,1:t_+1),fliplr(E_X(1,1:t_+1)-s_X(6,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')

% plot(t(1:t_+1),ReOrTrue(1,1:t_+1),'o--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','g','LineWidth',1)
% patch([t(1:t_+1),fliplr(t(1:t_+1))],[ReOrTrue(1,1:t_+1)+MVvTrue(1,1:t_+1),fliplr(ReOrTrue(1,1:t_+1)-MVvTrue(4,1:t_+1))],'g','FaceAlpha',0.2,'EdgeColor','none')

plot(t(1:t_+1),XTrue(1,1:t_+1),'--k','LineWidth',1)
grid on



xlabel('Time[Year]')
ylabel('Condition')
set(gca,'xtick',t)
xlim([0,T-1])
ylim([25,100])
xtickangle(45)
% h=legend('Inspection','$\pm2 \sigma_{Inspector}$','$\pm2 \sigma_{True_{Inspector}}$','Median','$\pm2 \sigma_{Model}$','$\pm1 \sigma_{Model}$','True Median','$\pm2 \sigma_{True}$','True State');
h=legend('Inspection','Inspection $- \mu_{Inspector}$', '$\pm2 \sigma_{Inspector}$','$\pm2 \sigma_{True_{Inspector}}$','Median','$\pm2 \sigma_{Model}$','$\pm1 \sigma_{Model}$','True State');

set(h,'Interpreter','latex')
set(h,'Position',[0.25 0.25 0.1 0.1])
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
ElementInd=1;
title(sprintf('Structure: %d, Poutre Element: %d',StructureInd,ElementInd))
hold off

figure(11)
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5],'PaperPositionMode','auto');
plot(t(1:t_+1),E_X(2,1:t_+1),'ro--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','r')
hold on
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(2,1:t_+1)+s_X(2,1:t_+1),fliplr(E_X(2,1:t_+1)-s_X(7,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(2,1:t_+1)+s_X(8,1:t_+1),fliplr(E_X(2,1:t_+1)-s_X(9,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')
% plot(t(1:t_+1),ExsmoothTrue(2,1:t_+1),'ro--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','g')
% patch([t(1:t_+1),fliplr(t(1:t_+1))],[ExsmoothTrue(2,1:t_+1)+MVvTrue(2,1:t_+1),fliplr(ExsmoothTrue(2,1:t_+1)-MVvTrue(5,1:t_+1))],'g','FaceAlpha',0.2,'EdgeColor','none')% just added
plot(t(1:t_+1),XbtrueTr(1,1:t_+1),'--k')
xlabel('Time')
ylabel('Speed')
set(gca,'xtick',t)
xlim([0,T-1])
xtickangle(45)
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
set(h,'Position',[0.25 0.25 0.1 0.1])
grid on
hold off

figure(12)
set(gcf,'Units','inches','Position',[0.5 0.5 5.5 3.5],'PaperPositionMode','auto');
plot(t(1:t_+1),E_X(3,1:t_+1),'ro--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','r')
hold on
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(3,1:t_+1)+2*s_X(3,1:t_+1),fliplr(E_X(3,1:t_+1)-2*s_X(3,1:t_+1))],'r','FaceAlpha',0.2,'EdgeColor','none')
patch([t(1:t_+1),fliplr(t(1:t_+1))],[E_X(3,1:t_+1)+s_X(3,1:t_+1),fliplr(E_X(3,1:t_+1)-s_X(3,1:t_+1))],'k','FaceAlpha',0.2,'EdgeColor','none')
% plot(t(1:t_+1),ExsmoothTrue(3,1:t_+1),'ro--','MarkerFaceColor',[1,1,1],'MarkerEdgeColor','g')
% patch([t(1:t_+1),fliplr(t(1:t_+1))],[ExsmoothTrue(3,1:t_+1)+MVvTrue(3,1:t_+1),fliplr(ExsmoothTrue(3,1:t_+1)-MVvTrue(3,1:t_+1))],'g','FaceAlpha',0.2,'EdgeColor','none')
plot(t(1:t_+1),XbbtrueTr(1,1:t_+1),'--k')
xlabel('Time')
ylabel('Acceleration')
set(gca,'xtick',t)
xlim([0,T-1])
xtickangle(45)
set(gca,'XTickLabel',num2str([StartDate:Step:EndDate]'))
grid on
hold off


end