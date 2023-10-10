function [Mtr,Mtrtr]=SpaceTransformationVec(Pn,Mv,MaxCond,MinCond)
% Space Transformation
% fct_ne=@(xts,nts) 1-((igamma(1/nts, xts.^nts))./(gamma(1/nts))); %unnormalized fct
% fct_n=@(xts,nts) gammainc(xts.^nts,1/nts);
invfct=@(yts,nts) (gammaincinv(abs(yts),1./nts)).^(1./nts);
% dfct_n=@(xts,nts)(nts*exp(-xts^nts))/gamma(1/nts);
n=2.^Pn;
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
L=length(Mv);
%ind_all = 1:L;
ind_max = find(Mv>MaxCond);
ind_upper = find(Mv>(MaxCond-MinCond)/2+MinCond);
ind_equal = find(Mv==(MaxCond-MinCond)/2+MinCond);
ind_min = find(Mv<MinCond);
ind_lower = find(Mv<(MaxCond-MinCond)/2+MinCond);
%a = setdiff(ind_all,[ind_max ind_upper ind_equal ind_min ind_lower]);
%ind_other = find(ismember(A,a));

Mtr = nan(1,L);
Mtrtr = nan(1,L);
if ~isempty(ind_max)
    ms=ones(1,length(ind_max));
    Mtr(ind_max)=invfct(ms,n);
    Mtrtr(ind_max)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(ind_max));
end
if ~isempty(ind_upper)
    ind_u0 = (Mv(ind_upper)==MaxCond);
    if ~isempty(ind_u0)
        Mv(ind_upper(ind_u0)) = MaxCond-0.0001;
    end
    ms=interp1([MinCond,MaxCond],[-1,1],Mv(ind_upper));
    ms((ms>0.995)) = 0.995;
    Mtr(ind_upper) = invfct(ms,n);
    ind_u = (Mtr(ind_upper)>MapAxes);
    if ~isempty(ind_u)
        Mtr(ind_upper(ind_u)) = MapAxes;
    end
    Mtrtr(ind_upper)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(ind_upper));
end
if ~isempty(ind_min)
    ms=-ones(1,length(ind_min));
    Mtr(ind_min)=-invfct(ms,n);
    ind_l = (Mtr(ind_min)<-MapAxes);
    if ~isempty(ind_l)
        Mtr(ind_min(ind_l)) = -MapAxes;
    end
    Mtrtr(ind_min)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(ind_min));
end
if ~isempty(ind_equal)
    Mtr(ind_equal)=0;
    Mtrtr(ind_equal)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(ind_equal));
end
if ~isempty(ind_lower)
    ind_l0 = (Mv(ind_lower)==MinCond);
    if ~isempty(ind_l0)
        Mv(ind_lower(ind_l0)) = MinCond+0.0001;
    end
    ms=interp1([MinCond,MaxCond],[-1,1],Mv(ind_lower));
    Mtr(ind_lower)=-invfct(ms,n);
    ind_il = (Mtr(ind_lower)<-MapAxes);
    if ~isempty(ind_il)
        Mtr(ind_lower(ind_il)) = -MapAxes;
    end
    Mtrtr(ind_lower)=interp1([-MapAxes,MapAxes],[25-MaxRange,100+MaxRange],Mtr(ind_lower));
end