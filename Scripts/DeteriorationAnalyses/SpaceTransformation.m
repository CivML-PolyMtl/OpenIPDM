function [Mtr,Mtrtr]=SpaceTransformation(Pn,Mv,MaxCond,MinCond)
% Space Transformation
% fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts);                    % Horizonntal axis [0,~]
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
MapAxes=invfct(0.999,n);
MaxRange=(MaxCond-MinCond)-(MaxCond-MinCond)/invfct(0.999,n);
L=length(Mv);
for i=1:L
    if Mv(i)>MaxCond
        ms=1;
        Mtr(i)=invfct(ms,n);
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[MinCond-MaxRange,MaxCond+MaxRange],Mtr(i));
    elseif Mv(i)>(MaxCond-MinCond)/2+MinCond
        if Mv(i)==MaxCond
            Mv(i)=MaxCond-0.0001;
        end
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        if ms>0.999
            ms=0.999;
        end
        Mtr(i)=invfct(ms,n);
        if Mtr(i)>MapAxes
            Mtr(i)=MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[MinCond-MaxRange,MaxCond+MaxRange],Mtr(i));
    elseif Mv(i)<MinCond
        ms=-1;
        Mtr(i)=-invfct(ms,n);
        if Mtr(i)<-MapAxes
            Mtr(i)=-MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[MinCond-MaxRange,MaxCond+MaxRange],Mtr(i));
    elseif Mv(i)==(MaxCond-MinCond)/2+MinCond
        Mtr(i)=0;
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[MinCond-MaxRange,MaxCond+MaxRange],Mtr(i));
    elseif Mv(i)<(MaxCond-MinCond)/2+MinCond
        if Mv(i)==MinCond
            Mv(i)=MinCond+0.0001;
        end
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        Mtr(i)=-invfct(ms,n);
        if Mtr(i)<-MapAxes
            Mtr(i)=-MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[MinCond-MaxRange,MaxCond+MaxRange],Mtr(i));
    else
        Mtr(i)=nan;
        Mtrtr(i)=nan;
    end
end