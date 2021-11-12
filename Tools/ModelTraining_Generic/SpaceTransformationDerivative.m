function dftr=SpaceTransformationDerivative(Pn,Mv,MaxCond,MinCond)
dev=@(yts,nts) (1./nts)*gampdf(abs(yts),1./nts,1)*(gammaincinv(abs(yts),1/nts))^((1/nts)-1);
n=2.^Pn;
L=length(Mv);
for i=1:L
    if Mv(i)>MaxCond
        ms=1;
        dftr(i)=dev(ms,n);
    elseif Mv(i)>(MaxCond-MinCond)/2+MinCond
        if Mv(i)==MaxCond
            Mv(i)=MaxCond-0.0001;
        end
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        if ms>0.995
            ms=0.995;
        end
        dftr(i)=dev(ms,1./n,1);
    elseif Mv(i)<MinCond
        ms=-1;
        dftr(i)=gampdf(abs(ms),1./n,1);
    elseif Mv(i)==(MaxCond-MinCond)/2+MinCond
        Mtr(i)=0;
        dftr(i)=1;
    elseif Mv(i)<(MaxCond-MinCond)/2+MinCond
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        Mtr(i)=-invfct(ms,n);
        dftr(i)=interp1([-2,2],[-12.5,137.5],Mtr(i));
    else
        Mtr(i)=nan;
        dftr(i)=nan;
    end
end
end