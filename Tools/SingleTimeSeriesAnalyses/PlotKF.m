% subplot(3,1,1)
figure(1)
x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(1,:)-permute(2*sqrt(Var(1,1,:)),[1,3,2]) fliplr(Ex(1,:)+permute(2*sqrt(Var(1,1,:)),[1,3,2]))];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Ex(1,:),'k','Linewidth',1)
    plot(timestamps,y,'o');
    errorbar(timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
elseif type==2 && originalspace
    y1=[ExKF(1,:)-Down2Std fliplr(ExKF(1,:)+Up2Std)];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,ExKF(1,:),'k','Linewidth',1)
    plot(timestamps,y,'o');  
    errorbar(timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
elseif type==2
    y1=[Exsmooth(1,:)-permute(2*sqrt(Vsmooth(1,1,:)),[1,3,2]) fliplr(Exsmooth(1,:)+permute(2*sqrt(Vsmooth(1,1,:)),[1,3,2]))];
    hold on
    patch(x1-1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Exsmooth(1,:),'k','Linewidth',1)
    plot(timestamps,y,'o');  
    errorbar(timestamps,y,Rlow',Rtop','LineStyle','none','CapSize',18,'Color','blue','Linewidth',1);
end
grid on
grid minor
box on
xlabel('Year','Interpreter', 'Latex')
ylabel('Condition','Interpreter', 'Latex')
ylim([25,100])

% subplot(3,1,2)
figure(2)
x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(2,:)-permute(2*sqrt(Var(2,2,:)),[1,3,2]) fliplr(Ex(2,:)+permute(2*sqrt(Var(2,2,:)),[1,3,2]))];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Ex(2,:),'k','Linewidth',1)
else
    y1=[Exsmooth(2,:)-permute(2*sqrt(Vsmooth(2,2,:)),[1,3,2]) fliplr(Exsmooth(2,:)+permute(2*sqrt(Vsmooth(2,2,:)),[1,3,2]))];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Exsmooth(2,:),'k','Linewidth',1)
end
grid on
grid minor
box on
xlabel('Year','Interpreter', 'Latex')
ylabel('Speed','Interpreter', 'Latex')

% subplot(3,1,3)
figure(3)
x1=[timestamps', fliplr(timestamps')]; 
if type==1
    y1=[Ex(3,:)-permute(2*sqrt(Var(3,3,:)),[1,3,2]) fliplr(Ex(3,:)+permute(2*sqrt(Var(3,3,:)),[1,3,2]))];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Ex(3,:),'k','Linewidth',1)
else
    y1=[Exsmooth(3,:)-permute(2*sqrt(Vsmooth(3,3,:)),[1,3,2]) fliplr(Exsmooth(3,:)+permute(2*sqrt(Vsmooth(3,3,:)),[1,3,2]))];
    hold on
    patch(x1,y1,1,'FaceColor','r','EdgeColor','none');
    alpha(.2); 
    plot(timestamps,Exsmooth(3,:),'k','Linewidth',1)
end
grid on
grid minor
box on
xlabel('Year','Interpreter', 'Latex')
ylabel('Acc','Interpreter', 'Latex')