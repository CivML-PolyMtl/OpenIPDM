function [sd, mu_v] = OnlineInference_One(Data_inspect, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM,init_estim, Ntr,MdataEngy)
%ONLINEINFERENCE_ONE Summary of this function goes here

%% lenghts
E = size(Data_inspect,1);
n   = size(Data_inspect,2);
n_x = size(init_x,1); %3;
I = size(UpdatedInspectorsData{1, 1}  ,1);



%lists for saving the evolution of the estimation

save_mu = zeros(2,I);
save_cov = zeros(5,5,I);
count = zeros(I,1);
mu_next = zeros(5,1);
cov_next = zeros(5,5);


%New process error
dt = 1;
Aki = A;
% Qki = PARAM(1).^2*[(dt^5)/20,(dt^4)/8,(dt^3)/6;(dt^4)/8,(dt^4)/4,(dt^3)/2;(dt^3)/6,(dt^3)/2,(dt^2)];
Qki = Q(PARAM(1));
At = blkdiag(Aki, eye(2));
Qt = blkdiag(Qki, zeros(2,2));

C = [1 0 0 1 1];


%% initialisation
mu_v0 = init_estim(1);
mu_sd_V0 = init_estim(2);

sd_v0 = init_estim(3);
sd_sd_v0 = init_estim(4);


%analyse
test = zeros(E,n);

n_e = 1;
LL = zeros(n_e,1);
save_epoch = zeros(n_e,I);


save_x_KF = cell(E,n);
save_V_KF = cell(E,n);


for i = 1:n_e
    
    save_x_KF = cell(E,n);
    save_V_KF = cell(E,n);
    
%     if firstrun == 1
        save_mu(1,:) = (UpdatedInspectorsData{1, 1}(:,2))'; %mu_v1
        save_cov(4,4,:) = mu_sd_V0^2;  %sig_v1
        
        save_mu(2,:) = (UpdatedInspectorsData{1, 1}(:,3).^2)' ; %mu_v2
        save_cov(5,5,:) = sd_sd_v0^2;  %sig_v2
%         
%     else 
%         save_mu(1,:) = mu_v0; %mu_v1
%         save_cov(4,4,:) = mu_sd_V0^2;  %sig_v1
%         
%         save_mu(2,:) = sd_v0^2 ; %mu_v2
%         save_cov(5,5,:) = sd_sd_v0^2;  %sig_v2
%     end
    %%Kalman filter
    for e = 1:E
        
        save_x_KF{e,1} = init_x(:,e);
        save_V_KF{e,1} = init_V(:,:,e);
        
        y =  Data_inspect(e,:);
        mu_now = cat(1,init_x(:,e),zeros(2,1));
        cov_now = blkdiag(init_V(:,:,e),zeros(2,2));
        
        for j = 1:n
            
            if isnan(y(j))
                %prediction
                cov_now = blkdiag(cov_now(1:n_x,1:n_x),zeros(2,2));
                mu_now = cat(1,mu_now(1:n_x),zeros(2,1));
                mu_next = At * mu_now;
                cov_next = At * cov_now * At' + Qt;
                cov_now = cov_next;
                mu_now = mu_next;
                
                %apply constraints on the deterioration speed
                if (mu_now(2)+2*sqrt(cov_now(2,2)))>0 && j~=1
                    Dc=[0 1 0  ;0 1 0];
                    d=[-50;0];
                    [mu_now(1:n_x),cov_now(1:n_x,1:n_x)]=KFConstraintsHandlingSeq(mu_now(1:n_x),cov_now(1:n_x,1:n_x),Dc,d,1);
                end
                
            else
                %index inspector
                num_ins   = Inspectorlabel(e,j);
                index_ins = find(UpdatedInspectorsData{1,1}==num_ins);
                count(index_ins) = 1 + count(index_ins);
                
                %prior
                cov_now = blkdiag(cov_now(1:n_x,1:n_x),zeros(2,2));
                mu_now = cat(1,mu_now(1:n_x),zeros(2,1));
                mu_now(4) = save_mu(1,index_ins);
                mu_now(5) = save_mu(2,index_ins);
                %%full covariance
                %                        cov_now(n_x+1:end,1:end) = save_cov(n_x+1:end,1:end,index_ins);
                %                        cov_now(1:n_x,n_x+1:end) = save_cov(1:n_x,n_x+1:end,index_ins);
                %%covariance mu/sig
                %                        cov_now(n_x+1:end,n_x+1:end) = save_cov(n_x+1:end,n_x+1:end,index_ins);
                %%covariance mu/x
                %                        cov_now(4,1:n_x) = save_cov(4,1:n_x,index_ins);
                %                        cov_now(1:n_x,4) = save_cov(1:n_x,4,index_ins);
                %%covariance sig/x
                %                        cov_now(5,1:n_x) = save_cov(5,1:n_x,index_ins);
                %                        cov_now(1:n_x,5) = save_cov(1:n_x,5,index_ins);
                %%sig_v et sig_v2hat
                cov_now(4,4) = save_cov(4,4,index_ins);
                cov_now(5,5) = save_cov(5,5,index_ins);
                
                %prediction
                mu_next = At * mu_now;
                cov_next = At * cov_now * At' + Qt;
                
                mu_next(5) = 0;
                cov_next(5,5) = save_mu(2,index_ins);
                
                %1st step update
                Gt = C*cov_next*C' ;
                Kt1 = cov_next*C'./Gt;
                yt = C*mu_next ;
                rt = y(j)- yt;
                cov_now = (eye(size(C,2)) - Kt1*C) * cov_next;
                mu_now  = mu_next + Kt1*rt;
                LL(i,1) =  LL(i,1) + log(normpdf(rt, 0, sqrt(Gt)));
                test(e,j) = yt;
                
                
                %2nd step update
                mu_vy           = mu_now(5); % mu_{t|t}^{V}
                sig_vy          = cov_now(5,5);
                
                myhat           = save_mu(2,index_ins);
                syhat           = save_cov(5,5,index_ins);
                
                mu_v2y          = mu_vy^2 + sig_vy; % mu_{t|t}^{V2}
                sig_v2y         = 2*sig_vy^2 + 4 * sig_vy * mu_vy^2;
                
                
                mu_v2           = myhat ;    % mu_{t|t-1}^{V2}
                cov_v2          = 3*syhat + 2*myhat^2;
                
                Kt2             = syhat/cov_v2;
                myhat           = myhat + Kt2*(mu_v2y - mu_v2); % mu_{t|t}^{V2hat}
                syhat           = syhat + (Kt2^2)*(sig_v2y - cov_v2);
                

                
                %apply constraints on the deterioration speed
                if sqrt(myhat)<1
                    Dc=[1 ; 1];
                    d=[0.5 ;10000];
                    [myhat,syhat]=KFConstraintsHandlingSeq(myhat,syhat,Dc,d,1);
                end
                
                %apply constraints on the deterioration speed
                if (mu_now(2)+2*sqrt(cov_now(2,2)))>0 && j~=1
                    Dc=[0 1 0 0 0 ;0 1 0 0 0];
                    d=[-50;0];
                    [mu_now,cov_now]=KFConstraintsHandlingSeq(mu_now,cov_now,Dc,d,1);
                end
                                
                %sauvegarde
                save_mu(1,index_ins) = mu_now(4);
                save_mu(2,index_ins) = myhat;
                save_cov(:,:,index_ins) = cov_now(:,:);
                save_cov(5,5,index_ins) = syhat;
                
            end
            
            save_x_KF{e,j} = mu_now;
            save_V_KF{e,j} = cov_now;
        
        end

        
    end
    save_epoch(i,:) = save_mu(2,:);
    
    for e = 1:E
        
        init_x(2,e) = 0;
        [~,MAxCondition]=SpaceTransformation(Ntr,100,100,25);
        num_ins   = Inspectorlabel(e,2);
        index_ins = find(UpdatedInspectorsData{1,1}==num_ins);
        NormObsValues = y(2) + save_mu(1,index_ins) ;
        NormObsValues(find(NormObsValues>MAxCondition))=MAxCondition;
        DifferenceObs=MAxCondition(end)-NormObsValues; 
        init_V(2,2,e) = (PARAM(4).^2).*(DifferenceObs)+(PARAM(6).^2);
    end
    
end

%Kalman smoother

save_x_KS = cell(E,n);
save_V_KS = cell(E,n);

for e = 1:E
    save_x_KS{e,end} = save_x_KF{e,end};
    save_V_KS{e,end} = save_V_KF{e,end};
end


for e = 1:E
    
    for j = 2:n

            cov_tt = save_V_KF{e,(n+1-j)}(1:3,1:3);
            mu_tt = save_x_KF{e,(n+1-j)}(1:3);
            
            mu_t1t = Aki*mu_tt;
            cov_t1t = Aki*cov_tt*Aki' + Qki;
            
            cov_t1T = save_V_KS{e,(n+2-j)}(1:3,1:3);
            mu_t1T = save_x_KS{e,(n+2-j)}(1:3);
            
            Jt = cov_tt*Aki'/cov_t1t;
            cov_tT = cov_tt + Jt*(cov_t1T-cov_t1t)*Jt';
            mu_tT = mu_tt + Jt*(mu_t1T - mu_t1t);
            
            if (mu_tT(2)+2*sqrt(cov_tT(2,2)))>0 && j~=1
                Dc=[0 1 0 ;0 1 0];
                d=[-50;0];
                [mu_tT,cov_tT]=KFConstraintsHandlingSeq(mu_tT,cov_tT,Dc,d,1);
            end
            
            save_x_KS{e,(n+1-j)} =  cat(1,mu_tT,zeros(2,1));
            save_V_KS{e,(n+1-j)} =  blkdiag(cov_tT,zeros(2,2));

        end
 
end



    for e = 1:E
        init_x(2,e) = 0;
        [~,MAxCondition]=SpaceTransformation(Ntr,100,100,25);
        num_ins   = Inspectorlabel(e,2);
        index_ins = find(UpdatedInspectorsData{1,1}==num_ins);
        NormObsValues = save_x_KS{e,2}(1) ;%+ save_mu(1,index_ins) ;
        NormObsValues(find(NormObsValues>MAxCondition))= MAxCondition;
        DifferenceObs = MAxCondition(end)-NormObsValues; % first time
        init_V(2,2,e) = (PARAM(4).^2).*(DifferenceObs)+(PARAM(6).^2);
    end
    
firstrun = 1;

[sd, mu_v] = OnlineInference_short(Data_inspect, Inspectorlabel, init_x, init_V, A, Q, UpdatedInspectorsData,PARAM, init_estim, MdataEngy);

end
  


