ExSmooth=zeros(size(x),'gpuArray');
VarSmooth=zeros(size(Var),'gpuArray');
if GPUCompute
    ExSmooth(:,:,TotalTimeSteps)=x(:,:,TotalTimeSteps);
    VarSmooth(:,:,:,TotalTimeSteps)=Var(:,:,:,TotalTimeSteps);
    for i=TotalTimeSteps-1:-1:1
        Xpred=pagefun(@mtimes,A,x(:,:,i));
        Vpred=pagefun(@mtimes,pagefun(@mtimes,A,Var(:,:,:,i)),pagefun(@transpose,A)) + Q;
        J=pagefun(@mtimes,pagefun(@mtimes,Var(:,:,:,i),pagefun(@transpose,A)),pagefun(@inv,Vpred));
        ExSmoothUpdate1=reshape(ExSmooth(:,:,i+1)-Xpred,3,1,size(ExSmooth(:,:,1),2));
        ExSmoothUpdate2=reshape(pagefun(@mtimes,J,ExSmoothUpdate1),3,size(ExSmooth(:,:,1),2));
        ExSmooth(:,:,i)=x(:,:,i)+ExSmoothUpdate2;
        VarSmooth(:,:,:,i)=Var(:,:,:,i)+pagefun(@mtimes,...
            pagefun(@mtimes,J,VarSmooth(:,:,:,i+1)-Vpred),pagefun(@transpose,J));
        ViIndex=find((2*sqrt(VarSmooth(2,2,:,i))+reshape(ExSmooth(2,:,i),1,1,...
            length(ExSmooth(2,:,i))))>GlobalCondData(2,4));
        if ~isempty(ViIndex) && GlobalCondData(2,1)
            xpredIN=ExSmooth(:,ViIndex,i);
            VpredIN=VarSmooth(:,:,ViIndex,i);
            d=[GlobalCondData(2,3);GlobalCondData(2,4)];
            D=[0 1 0;0 1 0];
            [xnewIN,VnewIN,Status]=KFConstraintsHandling(xpredIN,VpredIN,D,d,Nsigma);
            ExSmooth(:,ViIndex,i)=xnewIN;
            VarSmooth(:,:,ViIndex,i)=VnewIN;
        end
    end

else
    ExSmooth(:,TotalTimeSteps)=x(:,TotalTimeSteps);
    VarSmooth(:,:,TotalTimeSteps)=Var(:,:,TotalTimeSteps);
    for i=TotalTimeSteps-1:-1:1
        Xpred=A*x(:,i);
        Vpred=A*Var(:,:,i)*A' + Q;
        J=Var(:,:,i)*A'*inv(Vpred);
        ExSmooth(:,i)=x(:,i)+J*(ExSmooth(:,i+1)-Xpred);
        VarSmooth(:,:,i)=Var(:,:,i)+J*(VarSmooth(:,:,i+1)-Vpred)*J';
        ViIndex=find((2*sqrt(VarSmooth(2,2,i))+ExSmooth(2,i))>GlobalCondData(2,4));
        if ~isempty(ViIndex) && GlobalCondData(2,1)
            xpredIN=ExSmooth(:,i);
            VpredIN=VarSmooth(:,:,i);
            d=[GlobalCondData(2,3);GlobalCondData(2,4)];
            D=[0 1 0;0 1 0];
            [xnewIN,VnewIN,Status]=KFConstraintsHandling(xpredIN,VpredIN,D,d,Nsigma);
            ExSmooth(:,i)=xnewIN;
            VarSmooth(:,:,i)=VnewIN;
        end
    end
    
end