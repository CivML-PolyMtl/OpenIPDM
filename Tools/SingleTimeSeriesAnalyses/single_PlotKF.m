% subplot(3,1,1)

x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(1,:)-permute(2*sqrt(Var(1,1,:)),[1,3,2]) fliplr(Ex(1,:)+permute(2*sqrt(Var(1,1,:)),[1,3,2]))];
    hold(app.UIAxes,"on")
    patch(app.UIAxes,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes,timestamps,Ex(1,:),'k','Linewidth',1)
    plot(app.UIAxes,timestamps,y,'o');
    errorbar(app.UIAxes,timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
elseif type==2 && originalspace
    y1=[ExKF(1,:)-Down2Std fliplr(ExKF(1,:)+Up2Std)];
    hold(app.UIAxes,"on")
    patch(app.UIAxes,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes,timestamps,ExKF(1,:),'k','Linewidth',1)
    plot(app.UIAxes,timestamps,y,'o');  
    errorbar(app.UIAxes,timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
elseif type==2
    y1=[Exsmooth(1,:)-permute(2*sqrt(Vsmooth(1,1,:)),[1,3,2]) fliplr(Exsmooth(1,:)+permute(2*sqrt(Vsmooth(1,1,:)),[1,3,2]))];
    hold(app.UIAxes,"on")
    patch(app.UIAxes,x1-1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes,timestamps,Exsmooth(1,:),'k','Linewidth',1)
    plot(app.UIAxes,timestamps,y,'o');  
    errorbar(app.UIAxes,timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
end
grid(app.UIAxes,'on')
box(app.UIAxes, "on")
xlabel(app.UIAxes, 'Year','Interpreter', 'Latex')
ylabel(app.UIAxes, 'Condition','Interpreter', 'Latex')
ylim(app.UIAxes,[25,100])

% subplot(3,1,2)

x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(2,:)-permute(2*sqrt(Var(2,2,:)),[1,3,2]) fliplr(Ex(2,:)+permute(2*sqrt(Var(2,2,:)),[1,3,2]))];
    hold(app.UIAxes_2,"on")
    patch(app.UIAxes_2,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes_2,timestamps,Ex(2,:),'k','Linewidth',1)
else
    y1=[Exsmooth(2,:)-permute(2*sqrt(Vsmooth(2,2,:)),[1,3,2]) fliplr(Exsmooth(2,:)+permute(2*sqrt(Vsmooth(2,2,:)),[1,3,2]))];
    hold(app.UIAxes_2,"on")
    patch(app.UIAxes_2,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2); 
    plot(app.UIAxes_2,timestamps,Exsmooth(2,:),'k','Linewidth',1)
end
grid(app.UIAxes_2,'on')
box(app.UIAxes_2, "on")
xlabel(app.UIAxes_2, 'Year','Interpreter', 'Latex')
ylabel(app.UIAxes_2,'Speed','Interpreter', 'Latex')

% subplot(3,1,3)

x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(3,:)-permute(2*sqrt(Var(3,3,:)),[1,3,2]) fliplr(Ex(3,:)+permute(2*sqrt(Var(3,3,:)),[1,3,2]))];
    hold(app.UIAxes_3,"on")
    patch(app.UIAxes_3,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes_3,timestamps,Ex(3,:),'k','Linewidth',1)
else
    y1=[Exsmooth(3,:)-permute(2*sqrt(Vsmooth(3,3,:)),[1,3,2]) fliplr(Exsmooth(3,:)+permute(2*sqrt(Vsmooth(3,3,:)),[1,3,2]))];
    hold(app.UIAxes_3,"on")
    patch(app.UIAxes_3,x1,y1,1,'FaceColor','r','EdgeColor','none','FaceAlpha',0.2);
    plot(app.UIAxes_3,timestamps,Exsmooth(3,:),'k','Linewidth',1)
end
grid(app.UIAxes_3,'on')
box(app.UIAxes_3, "on")
xlabel(app.UIAxes_3, 'Year','Interpreter', 'Latex')
ylabel(app.UIAxes_3,'Acc','Interpreter', 'Latex')