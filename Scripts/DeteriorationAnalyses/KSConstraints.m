if InterventionVector(i+1)
    if ConstrainedKF && (ExSmooth(2,i)>ExSmooth(2,i+1) || ExSmooth(2,i)+2*sqrt(VarSmooth(2,2,i))>0)
        D=[0 1 0 0 0 0 ;0 1 0 0 0 0];
        d=[-75;ExSmooth(2,i+1)-0*sqrt(VarSmooth(2,2,i+1))];
        [ExSmooth(:,i),VarSmooth(:,:,i)]=KFConstraintsHandlingSeq(ExSmooth(:,i),VarSmooth(:,:,i),D,d,1);
    end
    if abs(ExSmooth(3,i))<abs(ExSmooth(3,i+1))
        D=[0 0 1 0 0 0 ;0 0 1 0 0 0 ];
        d=[-5;ExSmooth(3,i+1)+1*sqrt(VarSmooth(3,3,i));0;0;0;0];
        [ExSmooth(:,i),VarSmooth(:,:,i)]=KFConstraintsHandlingSeq(ExSmooth(:,i),VarSmooth(:,:,i),D,d,1);
    end
end
if ExSmooth(2,i)+2*sqrt(VarSmooth(2,2,i))>0 && ConstrainedKF && ~InterventionVector(i+1)
    d=[-50;0];
    D=[0 1 0 0 0 0;0 1 0 0 0 0];
    [ExSmooth(:,i),VarSmooth(:,:,i)]=KFConstraintsHandlingSeq(ExSmooth(:,i),VarSmooth(:,:,i),D,d,1);
end