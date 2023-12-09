function PlotTimeSeries(YearSteps,E_X,s_X,y,y_unbiased,Rtop,Rlow,AllInspData,InspectorLabel,InterventionVector,CondFigure,SpeedFigure,ColorCode, varargin)
args = varargin;
if ~isempty(args)
    InterventionCode = args{1};
end
if ColorCode==1
    CD1='ro--';
    CD2='r';
elseif ColorCode==2
    CD1='mo--';
    CD2='m';
elseif ColorCode==3
    CD1='bo--';
    CD2='b';
else
    CD1='ko--';
    CD2='k';
end
T=length(y);
IntTime=find(InterventionVector);
t=YearSteps(1):YearSteps(end);
if ~isempty(IntTime)
    t1=YearSteps(1):YearSteps(1)+(IntTime-2);
    t2=YearSteps(1)+IntTime-1:YearSteps(end);
else
    t1=t;
    t2=t;
end
cla(CondFigure);
% InspId=InspectorLabel(find(~isnan(y)));
for i=1:length(InspectorLabel)
    InspectorsID{i}=num2str(InspectorLabel(i));
end

plot(CondFigure,t1,E_X(1,1:length(t1)),CD1,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',CD2,'LineWidth',1)
hold(CondFigure);
patch(CondFigure,[t1,fliplr(t1)],[E_X(1,1:length(t1))+s_X(1,1:length(t1)),fliplr(E_X(1,1:length(t1))-s_X(4,1:length(t1)))],CD2,'FaceAlpha',0.2,'EdgeColor','none')
patch(CondFigure,[t1,fliplr(t1)],[E_X(1,1:length(t1))+s_X(5,1:length(t1)),fliplr(E_X(1,1:length(t1))-s_X(6,1:length(t1)))],'k','FaceAlpha',0.2,'EdgeColor','none')

if ~isempty(IntTime)
    plot(CondFigure,t2,E_X(1,length(t1)+1:end),CD1,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',CD2,'LineWidth',1)
    patch(CondFigure,[t2,fliplr(t2)],[E_X(1,length(t1)+1:end)+s_X(1,length(t1)+1:end),fliplr(E_X(1,length(t1)+1:end)-s_X(4,length(t1)+1:end))],CD2,'FaceAlpha',0.2,'EdgeColor','none')
    patch(CondFigure,[t2,fliplr(t2)],[E_X(1,length(t1)+1:end)+s_X(5,length(t1)+1:end),fliplr(E_X(1,length(t1)+1:end)-s_X(6,length(t1)+1:end))],'k','FaceAlpha',0.2,'EdgeColor','none')
    y_ver=25:100;%E_X(1,IntTime)-10:E_X(1,IntTime)+10;
    x_ver=(IntTime-1)*ones(1,length(y_ver))+YearSteps(1);
    plot(CondFigure,x_ver,y_ver,'k-.');
    y_ver2=25:100;%E_X(1,IntTime-1)-10:E_X(1,IntTime-1)+10;
    x_ver2=(IntTime-2)*ones(1,length(y_ver2))+YearSteps(1);
    plot(CondFigure,x_ver2,y_ver2,'k-.');
    for ip=1:length(IntTime)
        x_patch=[(IntTime(ip)-2)+YearSteps(1) (IntTime(ip)-2)+YearSteps(1) (IntTime(ip)-1)+YearSteps(1) (IntTime(ip)-1)+YearSteps(1)];
        y_patch=[min(y_ver) max(y_ver) max(y_ver) min(y_ver)];
        patch(CondFigure,x_patch,y_patch,'k','FaceAlpha',0.2,'EdgeColor','none')
        if InterventionCode == 1000 | InterventionCode == 2000 | InterventionCode == 3000
            text(CondFigure,(IntTime(ip)-1)+YearSteps(1), min(y_ver)+11, sprintf('Intervention \n Code: %d \n  (Inferred)',InterventionCode),'Rotation',00,'interpreter','latex','FontSize', 14);
        else
            text(CondFigure,(IntTime(ip)-1)+YearSteps(1), min(y_ver)+11, sprintf('Intervention \n Code: %d \n ',InterventionCode),'Rotation',00,'interpreter','latex','FontSize', 14);
        end
    end
end
plot(CondFigure,t,y,'o','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','b')
if ~isempty(AllInspData)
    plot(CondFigure,AllInspData(:,1),AllInspData(:,2),'*','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','k')
end
if ~isempty(y_unbiased)
    plot(CondFigure,t,y_unbiased,'*','MarkerFaceColor',[0,0,1],'MarkerEdgeColor','r')
    errorbar(CondFigure,t,y_unbiased,Rlow',Rtop','LineStyle','none','CapSize',20,'Color','blue','Linewidth',1)
else
    errorbar(CondFigure,t,y,Rlow',Rtop','LineStyle','none','CapSize',20,'Color','blue','Linewidth',1)
end
text(CondFigure,t(find(~isnan(y))), y(find(~isnan(y)))-Rlow(find(~isnan(y)))-1, InspectorsID);
set(CondFigure,'xtick',YearSteps)
xlim(CondFigure,[YearSteps(1),YearSteps(end)])
ylim(CondFigure,[25,100])
xtickangle(CondFigure,45)
set(CondFigure.YLabel,'String','Condition $\tilde{\mu}_{t|T}$','interpreter','latex')
hold(CondFigure,'off');


cla(SpeedFigure);
axis(SpeedFigure,'tight');
plot(SpeedFigure,t1,E_X(2,1:length(t1)),CD1,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',CD2);
hold(SpeedFigure);
patch(SpeedFigure,[t1,fliplr(t1)],[E_X(2,1:length(t1))+s_X(2,1:length(t1)),fliplr(E_X(2,1:length(t1))-s_X(7,1:length(t1)))],CD2,'FaceAlpha',0.2,'EdgeColor','none')
patch(SpeedFigure,[t1,fliplr(t1)],[E_X(2,1:length(t1))+s_X(8,1:length(t1)),fliplr(E_X(2,1:length(t1))-s_X(9,1:length(t1)))],'k','FaceAlpha',0.2,'EdgeColor','none')
if ~isempty(IntTime)
    plot(SpeedFigure,t2,E_X(2,length(t1)+1:T),CD1,'MarkerFaceColor',[1,1,1],'MarkerEdgeColor',CD2)
    patch(SpeedFigure,[t2,fliplr(t2)],[E_X(2,length(t1)+1:T)+s_X(2,length(t1)+1:T),fliplr(E_X(2,length(t1)+1:T)-s_X(7,length(t1)+1:T))],CD2,'FaceAlpha',0.2,'EdgeColor','none')
    patch(SpeedFigure,[t2,fliplr(t2)],[E_X(2,length(t1)+1:T)+s_X(8,length(t1)+1:T),fliplr(E_X(2,length(t1)+1:T)-s_X(9,length(t1)+1:T))],'k','FaceAlpha',0.2,'EdgeColor','none')
    y_ver=min(ylim(SpeedFigure)):0.001:0;(min(E_X(2,:))-max(2*s_X(8,:))):0.001:0;%E_X(1,IntTime)-10:E_X(1,IntTime)+10;
    x_ver=(IntTime-1)*ones(1,length(y_ver))+YearSteps(1);
    plot(SpeedFigure,x_ver,y_ver,'k-.');
    y_ver2=min(ylim(SpeedFigure)):0.001:0;%E_X(1,IntTime-1)-10:E_X(1,IntTime-1)+10;
    x_ver2=(IntTime-2)*ones(1,length(y_ver2))+YearSteps(1);
    plot(SpeedFigure,x_ver2,y_ver2,'k-.');
    for ip=1:length(IntTime)
        x_patch=[(IntTime(ip)-2)+YearSteps(1) (IntTime(ip)-2)+YearSteps(1) (IntTime(ip)-1)+YearSteps(1) (IntTime(ip)-1)+YearSteps(1)];
        y_patch=[min(y_ver) max(y_ver) max(y_ver) min(y_ver)];
        patch(SpeedFigure,x_patch,y_patch,'k','FaceAlpha',0.2,'EdgeColor','none')
    end
end
set(SpeedFigure,'xtick',YearSteps)
xlim(SpeedFigure,[YearSteps(1),YearSteps(end)])
xtickangle(SpeedFigure,45)
set(SpeedFigure.YLabel,'String','Speed $\tilde{\dot{\mu}}_{t|T}$','interpreter','latex')
hold(SpeedFigure,'off');
end

