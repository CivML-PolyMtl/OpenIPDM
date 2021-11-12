function [inspectors]=FIlteredInspectors(MdataEngy,inspectors,AW)
jj=1;
CID=sum(~cellfun(@isempty,MdataEngy),2);
    YS=NaN(1,sum(CID),AW);
    InspectorLabelS=NaN(1,sum(CID),AW);
    AllInspectors=cell(1,size(inspectors,1));
for i=1:length(CID)
    for j=1:CID(i)
        % stretch the y vector to 1 year basis
        while length(unique(MdataEngy{i,j}(:,3)))~=length(MdataEngy{i,j}(:,3))
            [~, ii] = unique(MdataEngy{i,j}(:,3), 'first');
            keep    = (diff([ii(:); length(MdataEngy{i,j}(:,3)) + 1]) > 1);
            Repindex  = ii(keep);
            if ~isempty(Repindex)
                RepInd=Repindex+1;
                MdataEngy{i,j}(RepInd,:)=[];
            end
        end
        NanDatabase=find(isnan(MdataEngy{i,j}(:,1)'));
        MdataEngy{i,j}(NanDatabase,:)=[];
        if ~isempty(MdataEngy{i,j})
            y=MdataEngy{i,j}(:,1)'; %y (database)
            if length(y)<3
                y=nan(5,1);
            else
                if max(abs(diff(y)))>15
                    y=nan(5,1);
                elseif length(diff(y))>2 && length(find(diff(y)>0))>=length(diff(y))-length(find(diff(y)>0)) && length(find(diff(y)>5))>=length(diff(y))-length(find(diff(y)>5))
                    y=nan(5,1);
                elseif length(diff(y))==2 && length(find(diff(y)>0))>=length(diff(y)) && length(find(diff(y)>5))>=length(diff(y))
                    y=nan(5,1);
                elseif MdataEngy{i,j}(1,3)<2000 
                    y=nan(5,1);
                else
                    % generate time sequence
                    yearly=MdataEngy{i,j}(1,3):MdataEngy{i,j}(end,3);
                    % find the years where that have observations
                    [~,~,iinsp]=intersect(MdataEngy{i,j}(:,3),yearly);
                    % stretch the observation vector over the time series
                    Diffy=nan(length(yearly),1);
                    Diffy(iinsp)=y;
                    y=Diffy';
                    % find inspectors for the time series
                    Insp=MdataEngy{i,j}(:,2)';
                    
                    InspectorLabel=nan(1,length(y));
                    InspectorLabel(iinsp)=Insp;
                end
            end
        else
            y=nan(5,1);
        end
        if sum(~isnan(y))>=3
            YS(1,jj,2:length(y)+1)=y;
            InspectorLabelS(1,jj,2:length(InspectorLabel)+1)=InspectorLabel;
            NaNIndex=find(~isnan(InspectorLabel));
            UILabel=unique(InspectorLabel(NaNIndex));
            for II=1:length(UILabel)
                ObsInspIdx=find(InspectorLabel==UILabel(II));
                InspectorIndex=find(inspectors(:,1)==UILabel(II));
                if isempty(AllInspectors{InspectorIndex})
                    AllInspectors{InspectorIndex}=zeros(sum(CID),AW);
                end
                AllInspectors{InspectorIndex}(jj,ObsInspIdx+1)=1;
            end
        end
        jj=jj+1;
    end
end

NaNTimeSeries=isnan(YS(1,:,2));
if ~isempty(NaNTimeSeries)
    YS(:,NaNTimeSeries,:)=[];
    for i=1:length(AllInspectors)
        if ~isempty(AllInspectors{i})
            AllInspectors{i}(NaNTimeSeries,:)=[];
        end
    end
end
ActiveInspectorsIndex=find(cellfun('isempty',AllInspectors));
inspectors(ActiveInspectorsIndex)=[];
end