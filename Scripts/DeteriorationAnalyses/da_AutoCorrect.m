function [YearTotal,y_Data,Re,InterventionVector,...
        InspectorsID,IntTime]=da_AutoCorrect(app,YearTotal,...
        y_Data,Re,InterventionVector,InspectorsID,InterventionDuration,...
        BudgetDuration,InterventionYear)
    DataInd=(~isnan(y_Data));
    y=y_Data(DataInd);
    Years=YearTotal(DataInd);
    R=Re(DataInd);
    
    % Initial check for interventions
    if  sum(InterventionVector)==0
        NoInt=1;
    else
        NoInt=0;
    end
   
    %% Outliers check
    % improve by 15 points
    if ~isempty(diff(y))
        if NoInt && max((diff(y)))>=15
            [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
                InterventionVector,InterventionDuration,BudgetDuration);
            if sum(InterventionVector)==0
                NoInt=1;
                while max((diff(y)))>=15
                    [y,R,InspectorsID,Years,YearTotal]=WeightedOutlier(y,R,InspectorsID,...
                        Years,YearTotal);
                end
            else
                NoInt=0;
            end
        end
    end
    % degrade by 15 points
    if ~isempty(diff(y))
        if NoInt && min((diff(y)))<=-15
            [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
                InterventionVector,InterventionDuration,BudgetDuration);
            if sum(InterventionVector)==0
                NoInt=1;
                while min((diff(y)))<=-15
                    [y,R,InspectorsID,Years,YearTotal]=WeightedOutlier(y,R,InspectorsID,...
                        Years,YearTotal);
                end
            else
                NoInt=0;
            end
        end
    end
    % staircase improvement by 5 points for half or more observations
    if ~isempty(diff(y))
    if NoInt && length(diff(y))>2 && 2*length(find(diff(y)>5))>=length(diff(y))
        [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
            InterventionVector,InterventionDuration,BudgetDuration);
        if sum(InterventionVector)==0
            NoInt=1;
            while length(diff(y))>2 && 2*length(find(diff(y)>5))>=length(diff(y))
                [y,R,InspectorsID,Years,YearTotal]=WeightedOutlier(y,R,InspectorsID,...
                    Years,YearTotal);
            end
        else
            NoInt=0;
        end
    end
    end
    
    
    % overall improvement by 15 points for time series with less than 8
    % Obbservations and sum of differences is negative
    if ~isempty(diff(y))
    if NoInt && max(y)-min(y)>15 && length(y)<8 && sum(diff(y))>-5
        [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
            InterventionVector,InterventionDuration,BudgetDuration);
        if sum(InterventionVector)==0
            NoInt=1;
            while max(y)-min(y)>15 && length(y)<8 && sum(diff(y))>-5
                [y,R,InspectorsID,Years,YearTotal]=WeightedOutlier(y,R,InspectorsID,...
                    Years,YearTotal);
            end
        else
            NoInt=0;
        end
    end
    end
    % overall improvement by 10 points for time series with less than 3
    % Obbservations
    if ~isempty(diff(y))
    if NoInt && max(y)-min(y)>10 && length(y)<3
        [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
            InterventionVector,InterventionDuration,BudgetDuration);
        if sum(InterventionVector)==0
            NoInt=1;
            while max(y)-min(y)>10 && length(y)<3
                [y,R,InspectorsID,Years,YearTotal]=WeightedOutlier(y,R,InspectorsID,...
                    Years,YearTotal);
            end
        else
            NoInt=0;
        end
    end
    end
    % check for an inntervention only
    if ~isempty(diff(y))
        if NoInt && max(diff(y))>=10
            [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
                InterventionVector,InterventionDuration,BudgetDuration);
        end
    end
     %% Outliers within Interventions check
     if ~NoInt && ~isempty(diff(y))
         y_parts=[];year_parts=[];R_parts=[];Insp_parts=[];YearTotal_parts=[];
         IndIntV=find(InterventionVector);
         [~,~,iyearv]=intersect(Years,YearTotal);
         for i=1:length(IndIntV)
             if i>1
                 y_before=[];
                 year_before=[];
                 R_before=[];
                 Insp_before=[];
                 YearTotal_before=[];
             else
                 IndBefore=~(IndIntV(i)<=iyearv);
                 y_before=y(IndBefore);
                 year_before=Years(IndBefore);
                 R_before=R(IndBefore);
                 Insp_before=InspectorsID(IndBefore);
                 YearTotal_before=YearTotal(1):YearTotal(IndIntV(i))-1;
                 while max(y_before)-min(y_before)>15
                     [y_before,R_before,Insp_before,year_before,YearTotal_before]=...
                         WeightedOutlier(y_before,R_before,Insp_before,year_before,...
                         YearTotal_before);
                 end
                 while min((diff(y_before)))<=-15
                     [y_before,R_before,Insp_before,year_before,YearTotal_before]=...
                         WeightedOutlier(y_before,R_before,Insp_before,year_before,...
                         YearTotal_before);
                 end
             end
             if i<length(IndIntV)
                 IndAfter=(~(IndIntV(i+1)<=iyearv)&IndIntV(i)<iyearv);
                 YearTotal_after=YearTotal(IndIntV(i)):YearTotal(IndIntV(i+1))-1;
             else
                 IndAfter=~(IndIntV(i)>iyearv);
                 YearTotal_after=YearTotal(IndIntV(i)):YearTotal(end);
             end
             y_after=y(IndAfter);
             year_after=Years(IndAfter);
             R_after=R(IndAfter);
             Insp_after=InspectorsID(IndAfter);
             
             
             
             while max(y_after)-min(y_after)>15
                 [y_after,R_after,Insp_after,year_after,YearTotal_after]=...
                     WeightedOutlier(y_after,R_after,Insp_after,year_after,...
                     YearTotal_after);
             end
             while min((diff(y_after)))<=-15
                 [y_after,R_after,Insp_after,year_after,YearTotal_after]=...
                     WeightedOutlier(y_after,R_after,Insp_after,year_after,...
                     YearTotal_after);
             end
             if length(IndIntV)>1
                 y_parts=[y_parts;y_before;y_after];
                 year_parts=[year_parts year_before year_after];
                 R_parts=[R_parts;R_before;R_after];
                 Insp_parts=[Insp_parts Insp_before Insp_after];
                 YearTotal_parts=[YearTotal_parts YearTotal_before YearTotal_after];
             end
             if ~isempty(InterventionYear)
                IntTime(i)=InterventionYear(i);
             end
         end
         if length(IndIntV)>1
             y=y_parts;
             Years=year_parts;
             R=R_parts;
             InspectorsID=Insp_parts;
             YearTotal=YearTotal_parts;
         else
             y=[y_before;y_after];
             Years=[year_before year_after];
             R=[R_before;R_after];
             InspectorsID=[Insp_before Insp_after];
             if isempty(YearTotal_before) || isempty(YearTotal_after)
                 YearTotal=[YearTotal_before YearTotal_after];
             else
                YearTotal=YearTotal_before(1):YearTotal_after(end);%[YearTotal_before YearTotal_after];
             end
         end
         
     else
         IntTime=0;
     end
    %% Time-series padding
    [~,~,iy]=intersect(Years,YearTotal);
    y_Data=nan(1,length(YearTotal));
    Re=y_Data;
    y_Data(iy)=y;
    Re(iy)=R;

end

function [InterventionVector,IntTime]=InterventionSearch(y,Years,YearTotal,...
        InterventionVector,InterventionDuration,BudgetDuration)
    if ~isempty(InterventionDuration)
        for i=1:size(InterventionDuration,1) % check if the year matches, and record the year
            AllYears=InterventionDuration(i,1):InterventionDuration(i,2);
            [~,Outlier]=max(diff(y));
            % time from the last observation before outlier : time of the outlier 
            OutlierYears=Years(Outlier):Years(Outlier+1);
            % earliest time after the last observation before outlier
            IntTime=min(intersect(AllYears,OutlierYears))+1;
            if ~isempty(IntTime)
                if IntTime>Years(Outlier+1)
                    IntTime=IntTime-1;
                end
                [~,IntInd]=intersect(IntTime,YearTotal);
                InterventionVector(IntInd)=1;
                break;
            end
        end
    else
        IntTime=[];
    end
    if sum(InterventionVector)==0 && ~isempty(BudgetDuration)
        for i=1:size(BudgetDuration,1) % check if the year matches, and record the year
            AllYears=BudgetDuration(i,1):BudgetDuration(i,2);
            [~,Outlier]=max(diff(y));
            % time from the last observation before outlier : time of the outlier 
            OutlierYears=Years(Outlier):Years(Outlier+1);
            % earliest time after the last observation before outlier
            IntTime=min(intersect(AllYears,OutlierYears))+1;
            if ~isempty(IntTime)
                if IntTime>Years(Outlier+1)
                    IntTime=IntTime-1;
                end
                [~,IntInd]=intersect(YearTotal,IntTime);
                InterventionVector(IntInd)=1;
                break;
            end
        end
    end
end

function [y,R,InspectorsID,Years,Duration]=WeightedOutlier(y,R,InspectorsID,Years,YearTotal)
WAvg=sum(y.*(1./R)./sum(1./R));
[~,Outlier]=max(abs(WAvg-y));
y(Outlier)=[];
R(Outlier)=[];
InspectorsID(Outlier)=[];
Years(Outlier)=[];
Duration=Years(1)-1:YearTotal(end);
end