function [xtb,Std,yOr,Rtop,Rlow,x_true]=BackTransformResults(y,Re,x,Std_Transformed,Pn,Xtrue,MaxCond,MinCond)
% Orgnize/Plot Data
ReOr=zeros(1,length(y));
invfct=@(yts,nts) (gammaincinv(abs(yts),1/nts))^(1/nts); % Horizonntal axis [0,~]
dfct_n=@(xts,nts)(nts*exp(-xts.^nts))/gamma(1/nts); 
n=2.^Pn;
MapAxes=invfct(0.999,n);
MaxExtendedAxes=(MaxCond-MinCond)-(MaxCond-MinCond)/invfct(0.999,n);


for i=1:length(x(1,:))
    xtb(1,i)=RevSpaceTransform(Pn,x(1,i),MaxCond,MinCond);                                  % Back-transform the condition (Model)
    yOr(i)=RevSpaceTransform(Pn,y(i),MaxCond,MinCond);                                      % Back-transform the condition (Observations)
    yTr(i)=interp1([MinCond-MaxExtendedAxes,MaxCond+MaxExtendedAxes],...
        [-MapAxes,MapAxes],y(i));                                           % Remap the observations from the transformation function values
    if x(1,i)>MaxCond+MaxExtendedAxes
        x(1,i)=MaxCond+MaxExtendedAxes;
    elseif x<MinCond-MaxExtendedAxes
        x(1,i)=MinCond-MaxExtendedAxes;
    end
    xTr(i)=interp1([MinCond-MaxExtendedAxes,MaxCond+MaxExtendedAxes],...
        [-MapAxes,MapAxes],x(1,i));                                         % Remap the model estimated condition from the transformation function values
    Std(1,i)=(RevSpaceTransform(Pn,2*Std_Transformed(1,i)+x(1,i),MaxCond,MinCond)...         % E[x(1)]+2*sigma - E[x(1)]
        -RevSpaceTransform(Pn,x(1,i),MaxCond,MinCond));
    Std(4,i)=(RevSpaceTransform(Pn,x(1,i),MaxCond,MinCond)-...                               % E[x(1)] - (E[x]-2*sigma)
        RevSpaceTransform(Pn,x(1,i)-2*Std_Transformed(1,i),MaxCond,MinCond));
    Std(5,i)=(RevSpaceTransform(Pn,1*Std_Transformed(1,i)+x(1,i),MaxCond,MinCond)-...        % E[x(1)]+1*sigma - E[x(1)]
        RevSpaceTransform(Pn,x(1,i),MaxCond,MinCond));
    Std(6,i)=(RevSpaceTransform(Pn,x(1,i),MaxCond,MinCond)-...                              % E[x(1)] - (E[x(1)]-1*sigma)
        RevSpaceTransform(Pn,x(1,i)-1*Std_Transformed(1,i),MaxCond,MinCond));
    ReOr(i)=(RevSpaceTransform(Pn,2*sqrt(Re(i))+y(i),MaxCond,MinCond)...                    % E[y]+2*sigma - (E[y]-2*sigma)
        -RevSpaceTransform(Pn,y(i)-2*sqrt(Re(i)),MaxCond,MinCond));
    Rtop(i)=RevSpaceTransform(Pn,2*sqrt(Re(i))+y(i),MaxCond,MinCond)...                     % E[y]+2*sigma - E[y]
        -RevSpaceTransform(Pn,y(i),MaxCond,MinCond);
    Rlow(i)=RevSpaceTransform(Pn,y(i),MaxCond,MinCond)-RevSpaceTransform(Pn,y(i)...         % E[y] - (E[y]-2*sigma)
        -2*sqrt(Re(i)),MaxCond,MinCond);
    xtb(2,i)=x(2,i)*dfct_n(xTr(i),n);                                    % E[x(2)] (original space)
    Std(2,i)=((2*Std_Transformed(2,i)+x(2,i)).*dfct_n(xTr(i),n))...                  % E[x(2)] + 2*sigma - E[x(2)]
        -(x(2,i).*dfct_n(xTr(i),n));
    Std(7,i)=((x(2,i).*dfct_n(xTr(i),n))-...                             % E[x(2)] - (E[x(2)] - 2*sigma)
        (x(2,i)-2*Std_Transformed(2,i)).*dfct_n(xTr(i),n));
    Std(8,i)=((1*Std_Transformed(2,i)+x(2,i)).*dfct_n(xTr(i),n)...             % E[x(2)] + 1*sigma - E[x(2)]
        -(x(2,i).*dfct_n(xTr(i),n)));
    Std(9,i)=((x(2,i).*dfct_n(xTr(i),n))-(x(2,i)...                      % E[x(2)] - (E[x(2)]-1*sigma)
        -1*Std_Transformed(2,i)).*dfct_n(xTr(i),n));
    % Remains without back transformation
    xtb(3,i)=x(3,i);%*(dfct_n(xTr(i),2^Pn))+x(2,i)^2*(d2fct_n(xTr(i),2^Pn));    % E[x(3)] (Transformed space)
    Std(3,i)=Std_Transformed(3,i);
    if ~isempty(Xtrue)
        x_true(1,i)=RevSpaceTransform(Pn,Xtrue(1,i),MaxCond,MinCond); 
        xTrTrue(i)=interp1([MinCond-MaxExtendedAxes,MaxCond+MaxExtendedAxes],[-MapAxes,MapAxes],Xtrue(1,i));  
        x_true(2,i)=Xtrue(2,i)*dfct_n(xTrTrue(i),n); 
        x_true(3,i)=Xtrue(3,i);
    else
        x_true=[];
    end
end
end