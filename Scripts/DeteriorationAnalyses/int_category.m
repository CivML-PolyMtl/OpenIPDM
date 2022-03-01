ElemAnalyses=0;
Ncurve=4;
d.Value = 0; 
YearsDuration=120;
YD=length(YearsDuration);
NumElms=length(Elms);
if NumElms~=0
    CondData=nan(NumElms,YD);
    SpeedData=nan(NumElms,YD);
    YearData=nan(NumElms,YD);
    QteData=nan(NumElms,YD);
    TypeData=nan(NumElms,YD);
    cecData=nan(NumElms,YD);
else
    CondData=nan(1,YD);
    SpeedData=nan(1,YD);
    YearData=nan(1,YD);
    QteData=nan(1,YD);
    TypeData=nan(1,YD);
    cecData=nan(1,YD);
end
for i=1:NumElms
    if ~NetCatAnalyses
        d=uiprogressdlg(app.MainWindow,'Title','Please Wait',...
            'Message',sprintf('Processing element %d/%d',i,length(Elms)));
        d.Value = (i/NumElms)/1.2;
        pause(0.001);
    end
    Elm = Elms(i);
    SlCatAnalyses=0;
    [~,~,~,~,~,IntervType]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses);
    if IntervType ~= 0
        SlCatAnalyses=1;
        [cond, speed, cec, InterventionYear, Qte, IntervType]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses);
    else
        cond=nan;
        speed=nan;
        InterventionYear=nan;
        Qte=nan;
        IntervType=nan;
        cec=nan;
    end
    CondData(i,1:length(cond))=cond;
    SpeedData(i,1:length(speed))=speed;
    YearData(i,1:length(InterventionYear))=InterventionYear;
    QteData(i,1:length(Qte))=Qte;
    TypeData(i,1:length(IntervType))=IntervType;
    cecData(i,1:length(cec))=cec;
end
if size(CondData,2)>1
    QteData = repmat(QteData,size(CondData,2),1);
    CondData = reshape(CondData,size(CondData,1)*size(CondData,2),1);
    SpeedData = reshape(SpeedData,size(SpeedData,1)*size(SpeedData,2),1);
    YearData = reshape(YearData,size(YearData,1)*size(YearData,2),1);
    TypeData = reshape(TypeData,size(TypeData,1)*size(TypeData,2),1);
    cecData = reshape(cecData,size(cecData,1)*size(cecData,2),1);
end
IntData = [CondData, SpeedData, cecData, YearData, QteData, TypeData];
ind_0 = (cecData==0);
IntData(ind_0,:)=[];

function [cond, speed, cec, InterventionYear, Qte, IntervType]=RunElementwiseSL(Elm,app,ElemAnalyses,CatAnalyses,StrucAnalyses,NetAnalyses,SlCatAnalyses)
    if SlCatAnalyses
        if app.InterventionCDFButton.Enable
            da_element();
            if find(InterventionVector)>1
                cond = xtb_v(1,find(InterventionVector)-1);
                speed = xtb_v(2,find(InterventionVector)-1);
                CEC_Data=nan(TimeWindow,1);
                CEC_Data(iaY)=CECData(UY);
                ind_cec = 1;
                int_vec_ind = find(InterventionVector);
                cec = CEC_Data(int_vec_ind-ind_cec);
                if length(IntervType)>length(int_vec_ind)
                    IntervType = IntervType(1:length(int_vec_ind));
                end
                for int_i = 1:length(int_vec_ind)
                    while isnan(CEC_Data(int_vec_ind(int_i)-ind_cec))
                        ind_cec = ind_cec + 1;
                        if int_vec_ind(int_i) - ind_cec <= 0
                            break;
                        else
                            cec(int_i) = CEC_Data(int_vec_ind(int_i)-ind_cec);
                        end
                    end
                end
            else
                cond=nan;
                speed=nan;
                InterventionYear=nan;
                IntervType=nan;
                cec=nan;
            end
        else
            cond=nan;
            speed=nan;
            InterventionYear=nan;
            IntervType = nan;
            Qte = nan;
            cec = nan;
        end
    else
        da_element();
        cond=nan;
        speed=nan;
        InterventionYear=nan;
        Qte = nan;
        cec = nan;
    end
end