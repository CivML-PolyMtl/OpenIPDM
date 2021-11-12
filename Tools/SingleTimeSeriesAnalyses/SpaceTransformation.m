function [Mtr,Mtrtr]=SpaceTransformation(Pn,Mv,MaxCond,MinCond)
% Space Transformation
% fct_ne=@(xts,nts) 1-((igamma(1/nts, xts.^nts))./(gamma(1/nts))); %unnormalized fct
% fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts);
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
if Pn==1
    MapAxes=2;
    MaxRange=1*75;
elseif Pn==2
    MapAxes=1.5;
    MaxRange=0.5*75;
elseif Pn==3
    MapAxes=1.2;
    MaxRange=0.2*75;
elseif Pn==4
    MapAxes=1.1;
    MaxRange=0.1*75;
elseif Pn==5
    MapAxes=1.02;
    MaxRange=0.02*75;
end
L=length(Mv);
for i=1:L
    if Mv(i)>MaxCond
        ms=1;
        Mtr(i)=invfct(ms,n);
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(i));
    elseif Mv(i)>(MaxCond-MinCond)/2+MinCond
        if Mv(i)==MaxCond
            Mv(i)=MaxCond-0.0001;
        end
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        if ms>0.995
            ms=0.995;
        end
        Mtr(i)=invfct(ms,n);
        if Mtr(i)>MapAxes
            Mtr(i)=MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(i));
    elseif Mv(i)<MinCond
        ms=-1;
        Mtr(i)=-invfct(ms,n);
        if Mtr(i)<-MapAxes
            Mtr(i)=-MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(i));
    elseif Mv(i)==(MaxCond-MinCond)/2+MinCond
        Mtr(i)=0;
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(i));
    elseif Mv(i)<(MaxCond-MinCond)/2+MinCond
        if Mv(i)==MinCond
            Mv(i)=MinCond+0.0001;
        end
        ms=interp1([MinCond,MaxCond],[-1,1],Mv(i));
        Mtr(i)=-invfct(ms,n);
        if Mtr(i)<-MapAxes
            Mtr(i)=-MapAxes;
        end
        Mtrtr(i)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(i));
    else
        Mtr(i)=nan;
        Mtrtr(i)=nan;
    end
end