function [Mtrv]=KernelSpaceTransform(Pn,Mv,MaxCond,MinCond)
% % Space Transformation
% fct_ne=@(xts,nts) 1-((igamma(1/nts, xts.^nts))./(gamma(1/nts))); %unnormalized fct
fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
% invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts);
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
L=length(Mv);
if Pn==1
    MapAxes=2;
    MaxExtendedAxes=(MaxCond-MinCond)/2;
elseif Pn==2
    MapAxes=1.5;
    MaxExtendedAxes=(MaxCond-MinCond)/4;
elseif Pn==3
    MapAxes=1.2;
    MaxExtendedAxes=(MaxCond-MinCond)/10;
elseif Pn==4
    MapAxes=1.1;
    MaxExtendedAxes=(MaxCond-MinCond)/20;
elseif Pn==5
    MapAxes=1.05;
    MaxExtendedAxes=(MaxCond-MinCond)/40;
else
    MapAxes=1;
    MaxExtendedAxes=0;
end
for i=1:L
    if Mv(i)>MaxCond+MaxExtendedAxes
        Mv(i)=MaxCond+MaxExtendedAxes;
    elseif Mv(i)<MinCond-MaxExtendedAxes
        Mv(i)=MinCond-MaxExtendedAxes;
    end
    Mv(i)=interp1([MinCond-MaxExtendedAxes,MaxCond+MaxExtendedAxes],[-MapAxes,MapAxes],Mv(i));
    if Mv(i)>MapAxes
        Mv(i)=MapAxes;
        Mtr=fct_n(Mv(i),n);
        Mtrv(i)=interp1([-1,1],[MinCond,MaxCond],Mtr);
    elseif Mv(i)>0
        Mtr=fct_n(Mv(i),n);
        Mtrv(i)=interp1([-1,1],[MinCond,MaxCond],Mtr);
    elseif Mv(i)<-MapAxes
        Mv(i)=-MapAxes;
        Mtr=-fct_n(Mv(i),n);
        Mtrv(i)=interp1([-1,1],[MinCond,MaxCond],Mtr);
    elseif Mv(i)==0
        Mtrv(i)=62.5;
    elseif Mv(i)<0
        Mtr=-fct_n(Mv(i),n);
        Mtrv(i)=interp1([-1,1],[MinCond,MaxCond],Mtr);
    elseif isnan(Mv(i))
        Mtrv(i)=nan;
    end
end
