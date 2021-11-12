function [Mtrv]=RevSpaceTransform(Pn,Mv)
% % Space Transformation
% fct_ne=@(xts,nts) 1-((igamma(1/nts, xts.^nts))./(gamma(1/nts))); %unnormalized fct
fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
% invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts);
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
L=length(Mv);
if Pn==1
    MapAxes=2;
    MaxExtendedAxes=37.5;
elseif Pn==2
    MapAxes=1.5;
    MaxExtendedAxes=18.75;
elseif Pn==3
    MapAxes=1.2;
    MaxExtendedAxes=7.5;
elseif Pn==4
    MapAxes=1.1;
    MaxExtendedAxes=3.75;
elseif Pn==5
    MapAxes=1.05;
    MaxExtendedAxes=1.875;
else
    MapAxes=1;
    MaxExtendedAxes=0;
end
for i=1:L
    if Mv(i)>100+MaxExtendedAxes
        Mv(i)=100+MaxExtendedAxes;
    elseif Mv(i)<25-MaxExtendedAxes
        Mv(i)=25-MaxExtendedAxes;
    end
    Mv(i)=interp1([25-MaxExtendedAxes,100+MaxExtendedAxes],[-MapAxes,MapAxes],Mv(i));
    if Mv(i)>MapAxes
        Mv(i)=MapAxes;
        Mtr=fct_n(Mv(i),n);
        Mtrv=interp1([-1,1],[25,100],Mtr);
    elseif Mv(i)>0
        Mtr=fct_n(Mv(i),n);
        Mtrv=interp1([-1,1],[25,100],Mtr);
    elseif Mv(i)<-MapAxes
        Mv(i)=-MapAxes;
        Mtr=-fct_n(Mv(i),n);
        Mtrv=interp1([-1,1],[25,100],Mtr);
    elseif Mv(i)==0
        Mtrv=62.5;
    elseif Mv(i)<0
        Mtr=-fct_n(Mv(i),n);
        Mtrv=interp1([-1,1],[25,100],Mtr);
    elseif isnan(Mv(i))
        Mtrv=nan;
    end
end
