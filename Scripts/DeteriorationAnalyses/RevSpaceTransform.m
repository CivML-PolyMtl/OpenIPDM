function [Mtrv]=RevSpaceTransform(Pn,Mv,MaxCond,MinCond)
% Space Transformation
fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);               % Vertical axis [0,1]   
invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts); % Horizonntal axis [0,~]
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
L=length(Mv);
MapAxes=invfct(0.999,n);
MaxExtendedAxes=(MaxCond-MinCond)-(MaxCond-MinCond)/invfct(0.999,n);

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
