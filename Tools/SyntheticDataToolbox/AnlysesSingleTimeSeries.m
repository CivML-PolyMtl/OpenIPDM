TimeSpanV=str2double(get(app.TimeSpan,'Value'));
TimeStepV=str2double(get(app.TimeStep,'Value'));
NumInspV=str2double(get(app.NumInsp,'Value'));
MinInspSTDV=str2double(get(app.MinInspSTD,'Value'));
MaxInspSTDV=str2double(get(app.MaxInspSTD,'Value'));
MaxCondV=str2double(get(app.MaxCondVal,'Value'));
MinCondV=str2double(get(app.MinCondVal,'Value'));
MaxBias=str2double(get(app.MaxBiasVal,'Value'));
MinBias=str2double(get(app.MinBiasVal,'Value'));
SigmaWV=str2double(get(app.SigmaW,'Value'));
Ncurve=str2double(get(app.TransN,'Value'));
Condition=get(app.CondSlider,'Value');
Speed=get(app.SpeedSlider,'Value');
Acc=get(app.AccSlider,'Value');
InterventionCheck=get(app.UnderwentInterventionsCheckBox,'Value');
InterventionsParam = app.InterventionsParam;
set(app.CondLbl, 'Text',num2str(Condition));
set(app.SpeedLbl, 'Text',num2str(Speed));
set(app.AccLbl, 'Text',num2str(Acc));
[x,Y,RejectTimeSeries,RejectReason]=GenerateSingleTimeSeries(Condition,...
    Speed,Acc,TimeSpanV,SigmaWV,Ncurve,MaxCondV,MinCondV,MaxInspSTDV,...
    MinInspSTDV,NumInspV,TimeStepV,InterventionCheck,InterventionsParam,...
    MaxBias,MinBias);
if RejectTimeSeries
    set(app.TSStatus, 'Text','Rejected');
    set(app.RJReason, 'Text',num2str(RejectReason));
else
    set(app.TSStatus, 'Text','Accepted');
    set(app.RJReason, 'Text','N.A.');
end
set(app.RJReason, 'Tooltip', ...
    ['[1]: Nan-Values in the Time-Series.', 10,...
    '[2]: Positive Speed Or High Acceleration.' 10,...
    '[3,4,5]: Exceeded allowed Speed.'  10,...
    '[6,7,8]: Exceeded allowed Acceleration.' 10,...
    '[9]: Init. Cond. is below the lowest allowed Init. Cond.' 10,...
    '[10, 11]: Exceeded allowed condition'])
TrueOriginal=RevSpaceTransform(Ncurve,x(1,:),MaxCondV,MinCondV);
figure(6)
plot(1964:(1963+length(TrueOriginal)),TrueOriginal);
title('True State Sample Data')
xlabel('Time')
ylabel('Structural Element Condition')
xlim([1964,1963+length(TrueOriginal)])
ylim([MinCondV,MaxCondV])
hold on
ObsOriginal=RevSpaceTransform(2,Y(2:end),MaxCondV,MinCondV);
plot(1965:(1964+length(ObsOriginal)),ObsOriginal,'*');
title('Observed State Sample Data')
xlabel('Time')
ylabel('Structural Element Condition')
xlim([1964,1963+length(TrueOriginal)])
ylim([MinCondV,MaxCondV])
grid minor
grid on
hold off