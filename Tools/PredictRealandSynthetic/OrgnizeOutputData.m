% Orgnize/Plot Data
ReOr=zeros(1,length(y));
for i=1:length(x(1,:))
    xtb(1,i)=VV_RevSpaceTransform(Pn,x(1,i));                                  % Back-transform the condition (Model)
    yOr(i)=VV_RevSpaceTransform(Pn,y(i));                                      % Back-transform the condition (Observations)
    yTr(i)=interp1([25-MaxExtendedAxes,100+MaxExtendedAxes],[-MapAxes,MapAxes],y(i));                              % Remap the observations from the transformation function values
    xTr(i)=interp1([25-MaxExtendedAxes,100+MaxExtendedAxes],[-MapAxes,MapAxes],x(1,i));                            % Remap the model estimated condition from the transformation function values
    MVv(1,i)=(VV_RevSpaceTransform(Pn,2*MVv(1,i)+x(1,i))...                    % E[x(1)]+2*sigma - E[x(1)]
        -VV_RevSpaceTransform(Pn,x(1,i)));
    MVv(4,i)=(VV_RevSpaceTransform(Pn,x(1,i))-...                              % E[x(1)] - (E[x]-2*sigma)
        VV_RevSpaceTransform(Pn,x(1,i)-2*s_Xsmooth(1,i)));
    MVv(5,i)=(VV_RevSpaceTransform(Pn,1*s_Xsmooth(1,i)+x(1,i))-...             % E[x(1)]+1*sigma - E[x(1)]
        VV_RevSpaceTransform(Pn,x(1,i)));
    MVv(6,i)=(VV_RevSpaceTransform(Pn,x(1,i))-...                              % E[x(1)] - (E[x(1)]-1*sigma)
        VV_RevSpaceTransform(Pn,x(1,i)-1*s_Xsmooth(1,i)));
    ReOr(i)=(VV_RevSpaceTransform(Pn,2*sqrt(Re(i))+y(i))...                    % E[y]+2*sigma - (E[y]-2*sigma)
        -VV_RevSpaceTransform(Pn,y(i)-2*sqrt(Re(i))));
    Rtop(i)=VV_RevSpaceTransform(Pn,2*sqrt(Re(i))+y(i))...                     % E[y]+2*sigma - E[y]
        -VV_RevSpaceTransform(Pn,y(i));
    Rlow(i)=VV_RevSpaceTransform(Pn,y(i))-VV_RevSpaceTransform(Pn,y(i)...         % E[y] - (E[y]-2*sigma)
        -2*sqrt(Re(i)));
    xtb(2,i)=x(2,i)*dfct_n(xTr(i),2^Pn);                                    % E[x(2)] (original space)
    MVv(2,i)=((2*MVv(2,i)+x(2,i)).*dfct_n(xTr(i),2^Pn))...                  % E[x(2)] + 2*sigma - E[x(2)]
        -(x(2,i).*dfct_n(xTr(i),2^Pn));
    MVv(7,i)=((x(2,i).*dfct_n(xTr(i),2^Pn))-...                             % E[x(2)] - (E[x(2)] - 2*sigma)
        (x(2,i)-2*s_Xsmooth(2,i)).*dfct_n(xTr(i),2^Pn));
    MVv(8,i)=((1*s_Xsmooth(2,i)+x(2,i)).*dfct_n(xTr(i),2^Pn)...             % E[x(2)] + 1*sigma - E[x(2)]
        -(x(2,i).*dfct_n(xTr(i),2^Pn)));
    MVv(9,i)=((x(2,i).*dfct_n(xTr(i),2^Pn))-(x(2,i)...                      % E[x(2)] - (E[x(2)]-1*sigma)
        -1*s_Xsmooth(2,i)).*dfct_n(xTr(i),2^Pn));
    Bias_org(i) = VV_RevSpaceTransform(Pn,y(i)-InpecBiase(i));                 % (E[y]-bias)
    xtb(3,i)=x(3,i);%*(dfct_n(xTr(i),2^Pn))+x(2,i)^2*(d2fct_n(xTr(i),2^Pn));    % E[x(3)] (Transformed space)
    if ~isempty(ReTrue)
        Xtrue(i)=VV_RevSpaceTransform(Pn,...                                   % The condition (Database) (Original Space)
            SynDatabaseState{ElementInd}(1,i+yearly(1)-1964));
        XtrueTr(i)=SynDatabaseState{ElementInd}(1,i+yearly(1)-1964);        % The condition (Database) (Transformed Space)
        xTrTrue(i)=interp1([25-MaxExtendedAxes,100+MaxExtendedAxes],[-MapAxes,MapAxes],ExsmoothTrue(1,i));
        xTrState(i)=interp1([25-MaxExtendedAxes,100+MaxExtendedAxes],[-MapAxes,MapAxes],XtrueTr(1,i));             % Remap the database condition from the transformation function values
        XbtrueTr(i)=SynDatabaseState{ElementInd}(2,i+yearly(1)-1964)...     % True speed (Original Space)
            .*dfct_n(xTrState(i),2^Pn);
        XbbtrueTr(i)=SynDatabaseState{ElementInd}(3,i+yearly(1)-1964);      % True Acceleration (Transformed space)
        Erbar(i)=abs(x(2,i)-SynDatabaseState{ElementInd}(2,i+yearly(1)-1964));    % abs( Model speed-True speed)
        Erbarbar(i)= abs(x(3,i)-XbbtrueTr(i));                              % abs(Model acc. - True acc.)
        AccVal(i)=abs(x(1,i)-XtrueTr(i));                                   % abs(Model Condition - Database Condition)
        Erbias(i)=(x(2,i)-SynDatabaseState{ElementInd}(2,i+yearly(1)-1964));    % (Model speed-True speed)
        Erbarbias(i)= (x(3,i)-XbbtrueTr(i));                              % (Model acc. - True acc.)
        AccValBias(i)=(x(1,i)-XtrueTr(i));                                   % (Model Condition - Database Condition)
        MVvTrue(1,i)=(VV_RevSpaceTransform(Pn,2*MVvTrue(1,i)+...               % E[x(1)] + 2*sigma - E[x(1)] (True parameters)
            ExsmoothTrue(1,i))-VV_RevSpaceTransform(Pn,ExsmoothTrue(1,i)));
        MVvTrue(4,i)=(VV_RevSpaceTransform(Pn,ExsmoothTrue(1,i))-...           % E[x(1)] - (E[x]-2*sigma) (True parameters)
            VV_RevSpaceTransform(Pn,ExsmoothTrue(1,i)-2*s_XsmoothTrue(1,i)));
        MVvTrue(2,i)=((2*MVvTrue(2,i)+ExsmoothTrue(2,i))...                 % E[x(2)] + 2*sigma - E[x(2)]
            .*dfct_n(xTrTrue(i),2^Pn))-(ExsmoothTrue(2,i).*dfct_n(xTrTrue(i),2^Pn));
        MVvTrue(5,i)=((ExsmoothTrue(2,i).*dfct_n(xTrTrue(i),2^Pn))-...      % E[x(2)] - (E[x(2)] - 2*sigma)
            (ExsmoothTrue(2,i)-2*s_Xsmooth(2,i)).*dfct_n(xTrTrue(i),2^Pn));
        AccSTD(i)=abs(AccVal(i))/(s_Xsmooth(1,i));                               % Normlized abs(Model Condition - Database Condition)/sigma
        XTrueParam(i)=VV_RevSpaceTransform(Pn,ExsmoothTrue(1,i));              % The condition (Model with True Params.) (Original Space)
        RtopTrue(i)=VV_RevSpaceTransform(Pn,2*sqrt(ReTrue(i))+y(i))...         % E[y] + 2*sigma - E[y]
            -VV_RevSpaceTransform(Pn,y(i));
        RlowTrue(i)=VV_RevSpaceTransform(Pn,y(i))-...                          % E[y] - (E[y]-2*sigma)
            VV_RevSpaceTransform(Pn,y(i)-2*sqrt(ReTrue(i)));
        if abs(x(1,i)-XtrueTr(i))/(2*s_Xsmooth(1,i))>=1
            PWithinCI(i)=0;
        else
            PWithinCI(i)=1;
        end
        if abs(ExsmoothTrue(1,i)-XtrueTr(i))/(2*s_XsmoothTrue(1,i))>=1
            PWithinCITrue(i)=0;
        else
            PWithinCITrue(i)=1;
        end
    else
        MAcc=[];
        AccVal=[];
        AccSTD=[];
        PWithinCI=[];
        PWithinCITrue=[];
        Erbar=[];
        Erbarbar=[];
        Erbias=[];
        Erbarbias=[];
        AccValBias=[];
        BoundViolation=[];
        XtrueValue=[];
        XValue=[];
        XbtrueTr=[];
        XbbtrueTr=[];
    end
end
y=yOr;
ReOr(find(isnan(ReOr)))=0;
Re=ReOr;
if ~isempty(ReTrue)
    MAcc=mean(AccVal);
    XValue=xtb;
    XtrueValue=Xtrue;
    BoundViolation=cumsum([0 diff(MVv(1,:)+xtb(1,:))])./(MVv(1,:)...
        +MVv(4,:));
    BoundViolation(BoundViolation<0)=0;
end