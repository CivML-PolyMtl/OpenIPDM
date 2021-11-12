function [AKr,X_ControlPoints]=KR_Prep(StrucAtt,KernelType,Kernel_l,Ncp)
AllAtt=reshape(gather(StrucAtt),size(StrucAtt,2),...
    size(StrucAtt,3));

Kr=1;X_ControlPoints=[];    
KernelType=string(KernelType);
for i=1:length(KernelType)
    % control points
    if strcmp(KernelType(i),'AAK')
        
        CatVec=min(AllAtt(:,i)):max(AllAtt(:,i));
        RepFactor=length(CatVec)/Ncp;
        if RepFactor>1
            StepSample=roundup(length(1:length(CatVec))/Ncp,1);
            IndexSample=1:StepSample:length(CatVec);
            while length(IndexSample)<Ncp
                LastCat= setdiff(CatVec,CatVec(IndexSample));
                NewInd=find(CatVec==LastCat(end));
                IndexSample=[IndexSample NewInd];
                IndexSample=sort(IndexSample);
            end
            FinalVec=CatVec(IndexSample);
        else
            RepFactor=ceil(1/RepFactor);
            repCatVec=repmat(CatVec,1,RepFactor);
            FinalVec=repCatVec(1:Ncp);
        end
        X_ControlPoints(:,i)=FinalVec';
    else
        X_ControlPoints(:,i)=linspace(min(AllAtt(:,i)),max(AllAtt(:,i)),Ncp);
    end
end
switch length(KernelType)
    case 2
        [Xat1,Xat2] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2));
         X_ControlPoints = [Xat1(:), Xat2(:)];
    case 3
        [Xat1,Xat2,Xat3] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2),...
            X_ControlPoints(:,3));
         X_ControlPoints = [Xat1(:), Xat2(:), Xat3(:)];
    case 4
        [Xat1,Xat2,Xat3,Xat4] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2),...
            X_ControlPoints(:,3),X_ControlPoints(:,4));
         X_ControlPoints = [Xat1(:), Xat2(:), Xat3(:), Xat4(:)];
    case 5
        [Xat1,Xat2,Xat3,Xat4,Xat5] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2),...
            X_ControlPoints(:,3),X_ControlPoints(:,4),X_ControlPoints(:,5));
         X_ControlPoints = [Xat1(:), Xat2(:), Xat3(:), Xat4(:), Xat5(:)];
    case 6
        [Xat1,Xat2,Xat3,Xat4,Xat5,Xat6] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2),...
            X_ControlPoints(:,3),X_ControlPoints(:,4),X_ControlPoints(:,5),...
            X_ControlPoints(:,6));
         X_ControlPoints = [Xat1(:), Xat2(:), Xat3(:), Xat4(:), Xat5(:), Xat6(:)];
    case 7
        [Xat1,Xat2,Xat3,Xat4,Xat5,Xat6,Xat7] = ndgrid(X_ControlPoints(:,1),X_ControlPoints(:,2),...
            X_ControlPoints(:,3),X_ControlPoints(:,4),X_ControlPoints(:,5),...
            X_ControlPoints(:,6),X_ControlPoints(:,7));
         X_ControlPoints = [Xat1(:), Xat2(:), Xat3(:), Xat4(:), Xat5(:), Xat6(:), Xat7(:)];
end
for i=1:length(KernelType)
    Krv=KernelFun(AllAtt(:,i),X_ControlPoints(:,i),Kernel_l(i),KernelType(i));
    Kr=Kr.*Krv;
end
% Kernel Function
AKr=Kr./sum(Kr,2);
end