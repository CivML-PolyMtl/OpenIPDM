function [Exsmooth,Vsmooth,s_Xsmooth,Status]=kalman_smoother(x,V,A,Q,Tdb,RunSmoother)

Exsmooth(:,Tdb) = x(:,Tdb);                                                 % expected value of the state at time (t=T)
Vsmooth(:,:,Tdb) = V(:,:,Tdb);                                              % variance of the state at time (t=T)
s_Xsmooth(:,Tdb) = sqrt(diag(V(:,:,Tdb)));                                  % std. of the state at time (t=T)
Status=[];                                                                  % Status of constraints
for ts=Tdb-1:-1:1
    xpred = A*x(:,ts);
    Vpred = A*V(:,:,ts)*A' + Q;                                             % Vpred = Cov[X(t+1) | t]
    J = V(:,:,ts) * A' * inv(Vpred);                                        % smoother gain matrix
    Exsmooth(:,ts) = x(:,ts) + J*(Exsmooth(:,ts+1) - xpred);
    Vsmooth(:,:,ts) = V(:,:,ts) + J*(Vsmooth(:,:,ts+1) - Vpred)*J';
    if (2*sqrt(Vsmooth(2,2,ts))+Exsmooth(2,ts))>0  && ts<=RunSmoother %&& ts~=1
        d=[-25;0];
        D=[0 1 0;0 1 0];
        [Exsmooth(:,ts),Vsmooth(:,:,ts),Status]=...
            KFConstraintsHandling(Exsmooth(:,ts),Vsmooth(:,:,ts),D,d,0);
    end
    if (-2*sqrt(Vsmooth(2,2,ts))+Exsmooth(2,ts))<-25 && ts<=RunSmoother % && ts~=1
        d=[-25;0];
        D=[0 1 0;0 1 0];
        [Exsmooth(:,ts),Vsmooth(:,:,ts),Status]=...
            KFConstraintsHandling(Exsmooth(:,ts),Vsmooth(:,:,ts),D,d,0);
    end
    s_Xsmooth(:,ts)=sqrt(diag(Vsmooth(:,:,ts)));
end
end