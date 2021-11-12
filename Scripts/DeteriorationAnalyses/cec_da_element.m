Elm = app.ElementListBox.Value;
ElmTypeInd=find(strcmp(app.BridgeData(:,end),Elm));
ElmType=app.BridgeData(ElmTypeInd(1),4);

%% Data
% Deterioration Info:
% 1: NOStruc % 2: DateInsp % 3: NoTrav % 4: Element % 5: NoElem % 6: Position
% 7: Materiau % 8: TypeElement % 9: A % 10: B % 11: C % 12: D % 13: CEC 
% 14: Inspection ID % 15: Element ID
VIColumns=[9,10,11,12]; CECCol=13; InspID=14; TimeCol=2;
DataInd=strcmp(Elm,app.BridgeData(:,15));
CECData=str2double(app.BridgeData(DataInd,CECCol));
InspData=str2double(app.BridgeData(DataInd,InspID));
TimeData=app.BridgeData(DataInd,TimeCol);
TimeData=datestr(datenum(TimeData));
YearData=str2num(TimeData(:,end-3:end));
TimeWindow=length(min(YearData):max(YearData))+1+round(app.ForecastYearsSlider.Value);
YearSteps=min(YearData)-1:max(YearData);
YearTotal=min(YearData)-1:max(YearData)+round(app.ForecastYearsSlider.Value);
[~,iaY]=intersect(YearSteps,YearData,'stable');
CEC_Data=nan(TimeWindow,1);
CEC_Data(iaY)=CECData;
cla(app.ElmCondition)
ylim(app.ElmCondition,[0.5,4.5])
plot(app.ElmCondition,YearTotal,CEC_Data,'ko','MarkerSize',12,'MarkerFaceColor', 'k')
