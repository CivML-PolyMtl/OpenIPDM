function [Mtrv]=RevSpaceTransformVec(Pn,Mv,MaxCond,MinCond)
% % Space Transformation
% fct_ne=@(xts,nts) 1-((igamma(1/nts, xts.^nts))./(gamma(1/nts))); %unnormalized fct
fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
% invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts);
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
L=length(Mv);
if Pn==1
    MapAxes=2;
    MaxRange=0.5*75;
elseif Pn==2
    MapAxes=1.5;
    MaxRange=0.25*75;
elseif Pn==3
    MapAxes=1.2;
    MaxRange=0.11*75;
elseif Pn==4
    MapAxes=1.1;
    MaxRange=0.05*75;
elseif Pn==5
    MapAxes=1.05;
    MaxRange=0.01*75;
else
    MapAxes=1;
    MaxRange=0*75;
end

Mtrv = nan(1,L);
ind_max = find(Mv>MaxCond+MaxRange);
ind_min = find(Mv<MinCond-MaxRange);

if ~isempty(ind_max)
    Mv(ind_max) = MaxCond+MaxRange;
end
if ~isempty(ind_min)
    Mv(ind_min) = MinCond-MaxRange;
end
Mv=interp1([MinCond-MaxRange,MaxCond+MaxRange],[-MapAxes,MapAxes],Mv );

ind_upper = find(Mv>MapAxes);
ind_zero_g = find(Mv>0);
ind_zero_e = find(Mv==0);
ind_zero_l = find(Mv<0);
ind_lower = find(Mv<-MapAxes);

if ~isempty(ind_upper)
    Mv(ind_upper)=MapAxes;
    Mtr=fct_n(Mv(ind_upper),n);
    Mtrv(ind_upper)=interp1([-1,1],[MinCond,MaxCond],Mtr);
end

if ~isempty(ind_zero_g)
    Mtr=fct_n(Mv(ind_zero_g),n);
    Mtrv(ind_zero_g)=interp1([-1,1],[MinCond,MaxCond],Mtr);
end

if ~isempty(ind_lower)
    Mv(ind_lower)=-MapAxes;
    Mtr=-fct_n(Mv(ind_lower),n);
    Mtrv(ind_lower)=interp1([-1,1],[MinCond,MaxCond],Mtr);
end

if ~isempty(ind_zero_e)
    Mtrv(ind_zero_e)=62.5;
end

if ~isempty(ind_zero_l)
    Mtr=-fct_n(Mv(ind_zero_l),n);
    Mtrv(ind_zero_l)=interp1([-1,1],[MinCond,MaxCond],Mtr);
end
end
