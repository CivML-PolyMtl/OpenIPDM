if InterventionVector(t+1)
    if ConstrainedKF && (ExpectedCond(2,t+1)+2*sqrt(Variance(2,2,t+1))>0 || ExpectedCond(2,t+1)<ExpectedCond(2,t))
        D=[0 1 0 0 0 0;0 1 0 0 0 0];
        d=[ExpectedCond(2,t)-0*sqrt(Variance(2,2,t));0];
        [ExpectedCond(:,t+1),Variance(:,:,t+1)]=KFConstraintsHandlingSeq(ExpectedCond(:,t+1),Variance(:,:,t+1),D,d,1);
    end
    if abs(ExpectedCond(3,t+1))>abs(ExpectedCond(3,t))
        D=[0 0 1 0 0 0;0 0 1 0 0 0];
        d=[1*sqrt(Variance(3,3,t));ExpectedCond(3,t)-1*sqrt(Variance(3,3,t))];
        [ExpectedCond(:,t+1),Variance(:,:,t+1)]=KFConstraintsHandlingSeq(ExpectedCond(:,t+1),Variance(:,:,t+1),D,d,1);
    end
end
if ConstrainedKF && (ExpectedCond(2,t+1)+2*sqrt(Variance(2,2,t+1)))>0 && ~InterventionVector(t+1) %&& t+1>=find(InterventionVector)
    D=[0 1 0 0 0 0;0 1 0 0 0 0];
    d=[-50;0];
    [ExpectedCond(:,t+1),Variance(:,:,t+1)]=KFConstraintsHandlingSeq(ExpectedCond(:,t+1),Variance(:,:,t+1),D,d,1);
end