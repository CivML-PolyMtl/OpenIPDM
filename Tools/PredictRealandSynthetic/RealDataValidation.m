function MdataEngyUpdate=RealDataValidation(MdataEngy,OldID,SSPDsorted,NewID)
idx = ~cellfun('isempty',MdataEngy);
OldStrIDs=unique(cellfun(@(c) c(1,4),MdataEngy(idx)));
idy = ~cellfun('isempty',SSPDsorted);
NewStrIDs=unique(cellfun(@(c) c(1,4),SSPDsorted(idy)));

OldStrCode=OldID(OldStrIDs,2);
NewStrCode=NewID(NewStrIDs,2);
% for i=1:size(MdataEngy,1)
%     OldStrIDs(i)=MdataEngy{i,1}(1,4);
% end
% for i=1:size(SSPDsorted,1)
%     NewStrIDs(i)=SSPDsorted{i,1}(1,4);
% end
MatchingRows2=[];


for i=1:length(OldStrIDs)
%     StrIDMap=find(NewStrIDs==OldStrIDs(i));
    StrIDMap=find(strcmp(NewStrCode,OldStrCode(i)));
    if ~isempty(StrIDMap)
        MatchingRows1(i,1)=StrIDMap;
        MdataEngyUpdate(i,:)=SSPDsorted(MatchingRows1(i),:);
    else
        MatchingRows2(i,1)=1;
    end
end
RemoveIndex1=find(MatchingRows1==0);
RemoveIndex2=find(MatchingRows2==1);
MdataEngyUpdate(RemoveIndex1,:)=[];
MdataEngy(RemoveIndex2,:)=[];
OldLen = cellfun(@(x) size(x,1), MdataEngy(:,1));
NewLen = cellfun(@(x) size(x,1), MdataEngyUpdate(:,1));
Delta=NewLen-OldLen;
NoUpdateInd=find(Delta==0 | Delta>1);
MdataEngyUpdate(NoUpdateInd,:)=[];
end