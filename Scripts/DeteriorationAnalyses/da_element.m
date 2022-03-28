if app.IntServiceLifeButton.Enable || SlCatAnalyses 
    NumYearsSlider=120;
    ColorCode=1;
elseif ElemAnalyses
    NumYearsSlider=round(app.ForecastYearsSlider.Value);
    ColorCode=1;
elseif CatAnalyses || CatReport
    NumYearsSlider=round(app.ForecastYearsSlider_2.Value);
elseif StrucAnalyses
    NumYearsSlider=round(app.ForecastYearsSliderStruc.Value);
elseif NetAnalyses
    NumYearsSlider=round(app.ForecastYearsSliderNet.Value);
end

ElmTypeInd=find(strcmp(app.BridgeData(:,end),Elm));
ElmType=app.BridgeData(ElmTypeInd(1),4);

%% Model Parameters
% Deteriorration Model
FN=struct2cell(app.ModelParam);
ModelParam=FN{1};
NAParam=FN{2};
NAIntParam=FN{3};
ETInd=find(strcmp(ModelParam(:,1),ElmType));
PriorParam=ModelParam{ETInd,2};
if ~isempty(PriorParam)
    InspectorsParam=ModelParam{ETInd,3};
    RegressionModel=ModelParam{ETInd,4};
    SA_Index=ModelParam{ETInd,5};
    MaterialIndex=ModelParam{ETInd,9}.Material;
    ElementTypeIndex=ModelParam{ETInd,9}.TypeElement;
    
    Ncurve = ModelParam{ETInd,10};
    app.curve_param = Ncurve;
    % Intervention Parameters
    IntParam=ModelParam{ETInd,8};
    ExIntParam=ModelParam{ETInd,6};
    VarIntParam=ModelParam{ETInd,7};
    
    %% Data
    % Deterioration Info:
    % 1: NOStruc % 2: DateInsp % 3: NoTrav % 4: Element % 5: NoElem % 6: Position
    % 7: Materiau % 8: TypeElement % 9: A % 10: B % 11: C % 12: D % 13: CEC
    % 14: Inspection ID % 15: Element ID
    VIColumns=[9,10,11,12]; CECCol=13; InspID=14; TimeCol=2;
    DataInd=strcmp(Elm,app.BridgeData(:,15));
    VIData=str2double(app.BridgeData(DataInd,VIColumns));
    y_values=VIData(:,1)+0.75*VIData(:,2)+0.5*VIData(:,3)+0.25*VIData(:,4);
    if sum(~isnan(y_values))>0
    [~,y_values_trans]=SpaceTransformation(Ncurve,y_values,100,25);
    CECData=str2double(app.BridgeData(DataInd,CECCol));
    InspData=str2double(app.BridgeData(DataInd,InspID));
    TimeData=app.BridgeData(DataInd,TimeCol);
    IndNanV=isnan(y_values);
    if sum(IndNanV)>0
        y_values_trans(IndNanV)=[];
        CECData(IndNanV)=[];
        InspData(IndNanV)=[];
        TimeData(IndNanV)=[];
    end
    TimeData=datestr(datenum(TimeData),'mm/dd/yyyy');
    [YearData,UY]=unique(str2num(TimeData(:,end-3:end)));
    TimeWindow=length(min(YearData):max(YearData))+1+NumYearsSlider;
    YearSteps=min(YearData)-1:max(YearData);
    YearTotal=min(YearData)-1:max(YearData)+NumYearsSlider;
    if ~ElemAnalyses
        AllBridgeDates=(datestr(datenum(app.BridgeData(:,2))));
        AllBridgeYears=unique(str2num(AllBridgeDates(:,end-3:end)));
        if max(YearData)<max(AllBridgeYears)
            YearTotal=min(YearData)-1:max(AllBridgeYears)+NumYearsSlider;
            TimeWindow=length(YearTotal);
        end
    end
    [~,iaY]=intersect(YearSteps,YearData,'stable');
    y_Data=nan(TimeWindow,1);
    y_Data(iaY)=y_values_trans(UY);
    InspData=InspData(UY);
    
    % Inspectors
    % Database Columns
    % 1: Inspection ID % 2: Inspector ID
    AllInspections=app.InspectorsData;
    i_ind=1;
    for i=1:length(InspData)
        if ~isempty(find(AllInspections(:,1)==InspData(i)))
            ID_Insp(i_ind)=find(AllInspections(:,1)==InspData(i));
            i_ind=i_ind+1;
        else
            y_Data(iaY(i))=nan;
            iaY(i)=[];
        end
    end
    Inspectors=AllInspections(ID_Insp,2);
    
    if ElemAnalyses
        DescriptionInd=find(DataInd);
        app.TypeElemLabel.Text=app.BridgeData(DescriptionInd(1),8);
        app.SpanNoElemLabel.Text=app.BridgeData(DescriptionInd(1),3);
        app.PositionElemLabel.Text=app.BridgeData(DescriptionInd(1),6);
        app.MaterialElemLabel.Text=app.BridgeData(DescriptionInd(1),7);
    end
    
    % Structural Attributes
    % Database Columns
    % 1: Bridge ID % 2: Structure Type % 3: Munipalite % 4: LAT % 5: Long
    % 6: Longtotale % 7: LongTablier % 8: LargHorstout % 9: LargCarrossable
    % 10: SuperfTablier % 11: AnConst % 12: Statut % 13: DJMA % 14: x_camions
    % 15: NBREVoies.
    StrucAtt=app.BridgeAttributes;
    Age=YearData(1)-str2double(StrucAtt(11));
    
    SelectedCat=app.BridgeData(DataInd,7);
    IndexConv=find(strcmp(MaterialIndex,SelectedCat(1)));
    if length(MaterialIndex)<3 && length(ElementTypeIndex)>length(MaterialIndex)
        SelectedCat=app.BridgeData(DataInd,8);
        IndexConv=find(strcmp(ElementTypeIndex,SelectedCat(1)));
    end
    if isempty(IndexConv)
        IndexConv=1;
    end
    
    SA1=str2double(StrucAtt(:,4:10));
    SA2=str2double(StrucAtt(:,13:15));
    StructuralAttributes=[IndexConv Age SA1 SA2];
    
    % Element Details
    % Database Columns
    % 1: Bridge ID % 2: NoTrav % 3: Element ID % 4: Element Category % 5: No. Element
    % 6: Position % 7: Materiau % 8: TypeElement % 9: Importance % 10: QTE Calculee
    ID_Qte=strcmp(Elm,app.Bridge_Details(:,3));
    Qte=str2double(app.Bridge_Details(ID_Qte,10));
    if isnan(Qte)
        QT_Ind=strcmp(app.Bridge_Details(:,4),ElmType);
        QT_Vals=app.Bridge_Details(QT_Ind,10);
        Qte=nanmean(str2double(QT_Vals));
        if isnan(Qte)
            Qte=1;
        end
    end
    Importance=app.Bridge_Details(ID_Qte,9);
    if strcmp(Importance,'Secondaire')
        ImportanceInd=2;
    else
        ImportanceInd=1;
    end
    
%     % Inspectors
%     % Database Columns
%     % 1: Inspection ID % 2: Inspector ID
%     AllInspections=app.InspectorsData;
%     i_ind=1;
%     for i=1:length(InspData)
%         if ~isempty(find(AllInspections(:,1)==InspData(i)))
%             ID_Insp(i_ind)=find(AllInspections(:,1)==InspData(i));
%             i_ind=i_ind+1;
%         end
%     end
%     Inspectors=AllInspections(ID_Insp,2);
    
    % Interventions file
    % Database Columns
    % 1: Bridge ID % 2: Intervention Start Year %3: Intervention End Year 
    % 4: Intervention Code % 5: Element Category % 6: Element ID
    ID_Int=strcmp(Elm,app.Bridge_Int(:,6));
    InterventionType=str2double(app.Bridge_Int(ID_Int,4));
    InterventionYear=unique(str2double(app.Bridge_Int(ID_Int,2)));
    InterventionDuration=unique(str2double(app.Bridge_Int(ID_Int,2:3)),'rows');
    if isempty(InterventionType)
        InterventionType=0;
    end
    
    % Budget file
    % 1: Bridge ID % 2: Start Year % 3: End Year
    BudgetDuration=str2double(app.Bridge_Cost(:,2:3));
    if isempty(BudgetDuration)
        InterventionsCost=0;
    end
    %% SSM | SSM-KR Configuration
    ModelConfig();
    %% Run Analyses
    % SSM | SSM-KR
    if sum(~isnan(y_Data))>0
        if FindType
            find_InterventionType();
        else
            [Ex, VarKF, loglikelihood, Exsmooth, Vsmooth, InterventionMu_Network,...
                InterventionVar_Network,~,~,~,~,LifeSpanInt]=KF_KS(y_Data,A,C,Q,Q_r,...
                Re,Be,PriorParam,init_x,init_V,Ncurve,ConstrainedKF,InterventionCheck,...
                InterventionVector,InterventionMu_Network,InterventionVar_Network);
        end
        for j=1:size(Vsmooth,3)
            Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
        end
        [xtb_v,Std_v,~,~,~,~]=BackTransformResults(y_Data,Re,Exsmooth,Std,Ncurve,y_Data-reshape(Be,1,[]),100,25);
        if ElemAnalyses
            if app.IntServiceLifeButton.Enable==1
                cla(app.ElmSpeed);
                xticks(app.ElmSpeed,'auto');
                LifeSpanInt=diff(LifeSpanInt);
                TotLifSpan=LifeSpanInt(LifeSpanInt>=10^-4);
                IndServiceLife=(LifeSpanInt>=10^-4);
                x_SL=InterventionYear:InterventionYear+length(LifeSpanInt)-1;
                x_SL=x_SL(IndServiceLife);
                if ~isempty(x_SL)
                    xlim(app.ElmSpeed,[x_SL(1) x_SL(end)]);
                    plot(app.ElmSpeed,x_SL,TotLifSpan);
                    [~,IndMaxSL]=max(TotLifSpan);
                    hold(app.ElmSpeed)
                    plot(app.ElmSpeed,[x_SL(IndMaxSL) x_SL(IndMaxSL)],[min(ylim(app.ElmSpeed)) max(ylim(app.ElmSpeed))],'-.');
                    text(app.ElmSpeed,x_SL(IndMaxSL)-1, max(ylim(app.ElmSpeed))-0.9*max(ylim(app.ElmSpeed)), sprintf('%d',x_SL(IndMaxSL)),'Rotation',90);
                    ylabel(app.ElmSpeed,'Pr. of Return (Condition)'); 
                    hold(app.ElmSpeed)
                else
                    msgbox('Analysis is not feasible for this case.', 'Analyses can not be performed','warn');
                end
            else
                ylabel(app.ElmSpeed,'Speed');
%                 for j=1:size(Vsmooth,3)
%                     Std(:,j)=sqrt(diag(Vsmooth(:,:,j)));
%                 end
                if ~isempty(y_Data_before)
                    y_Data_before(:,2)=RevSpaceTransform(Ncurve,y_Data_before(:,2),100,25);
                end
                [xtb,Std,yOr,Rtop,Rlow,yOr_unbiased]=BackTransformResults(y_Data,Re,Exsmooth,Std,Ncurve,y_Data-reshape(Be,1,[]),100,25);
                PlotTimeSeries(YearTotal,xtb,Std,yOr,yOr_unbiased,Rtop,Rlow,y_Data_before,InspectorIDLabel_y,InterventionVector,app.ElmCondition,app.ElmSpeed,ColorCode);
            end
            if sum(InterventionVector)>0
                app.IntServiceLifeButton.Enable=1;
            else
                app.IntServiceLifeButton.Enable=0;
            end
        end
    else
        [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Be,Qte,MisStart]=NoAnalyses(app);
        IntervType=0;
        if ElemAnalyses
            msgbox('The selected structural element has no inspection data, or it has no model associated with it.', 'Analyses can not be performed','warn');
        end
    end
    else
        [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Be,Qte,MisStart]=NoAnalyses(app);
        IntervType=0;
        if ElemAnalyses
            msgbox('The selected structural element has no inspection data, or it has no model associated with it.', 'Analyses can not be performed','warn');
        end
    end
    
else
    [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Be,Qte,MisStart]=NoAnalyses(app);
    IntervType=0;
    if ElemAnalyses
        msgbox('The selected structural element has no inspection data, or it has no model associated with it.', 'Analyses can not be performed','warn');
        cla(app.ElmSpeed);
        cla(app.ElmCondition);
    end
end
clc

function [Exsmooth,Vsmooth,YearTotal,y_Data,Re,Be,Qte,MisStart]=NoAnalyses(app)
    YearsDuration=app.TotalYearsDuration;
    NumNanYears=length(YearsDuration);
    Exsmooth=nan(3,NumNanYears);
    Vsmooth=nan(3,3,NumNanYears);
    YearTotal=YearsDuration;
    y_Data=nan(1,NumNanYears);
    Re=nan(1,NumNanYears);
    Be=nan(1,NumNanYears);
    Qte=0;
    MisStart=YearTotal(1);
end