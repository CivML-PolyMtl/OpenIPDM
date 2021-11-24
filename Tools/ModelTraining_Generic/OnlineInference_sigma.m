function [sd] = OnlineInference_sigma(Data_inspect, Inspectorlabel, init_x, init_V, A,Q, UpdatedInspectorsData,PARAM, init_estim, MdataEngy);

%ONLINEINFERENCE_SHORT Summary of this function goes here
%% version 4 sigma only (size hidden state vector 4)

%lenghts
E = size(Data_inspect,1);
n   = size(Data_inspect,2);
n_x = size(init_x,1);
I = size(UpdatedInspectorsData{1, 1}  ,1);


%lists for saving the evolution of the estimation
save_var_mu = zeros(E,I);
save_var_sig = zeros(E,I);
save_mu = zeros(1,I);
save_cov = zeros(4,4,I);

count = zeros(I,1);

mu_next = zeros(4,1);
cov_next = zeros(4,4);


%New process error
dt = 1;
Aki = A;
Qki = Q(PARAM(1));

At = blkdiag(Aki, 1);
Qt = blkdiag(Qki, 0);

C = [1 0 0 1];

%initial estimate

sd_v0 = init_estim(1);
sd_sd_v0 = init_estim(2);

save_mu(1,:) = (UpdatedInspectorsData{1, 1}(:,end).^2)' ;
save_cov(4,4,:) = sd_sd_v0^2;


test = zeros(E,n);

n_e = 1;
LL = zeros(n_e,1);
save_epoch = zeros(n_e,I);

for i = 1:n_e
    %           %%Kalman filter
    for e = 1:E
        
        
        y =  Data_inspect(e,:);
        % Kr initial vector
        mu_now = cat(1,init_x(:,e),0);
        cov_now = blkdiag(init_V(:,:,e),0);
        
        for j = 1:n
            
            if isnan(y(j))
                %prediction
                cov_now = blkdiag(cov_now(1:n_x,1:n_x),0);
                mu_now = cat(1,mu_now(1:n_x),0);
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
                cov_now = blkdiag(cov_now(1:n_x,1:n_x),0);
                mu_now = cat(1,mu_now(1:n_x),0);
                mu_now(4) = save_mu(1,index_ins);
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
                
                
                %prediction
                mu_next = At * mu_now;
                cov_next = At * cov_now * At' + Qt;
                
                mu_next(4) = 0;
                cov_next(4,4) = save_mu(1,index_ins);
                
                %1st step update
                Gt = C*cov_next*C' ;
                Kt1 = cov_next*C'/Gt;
                yt = C*mu_next;
                rt = y(j)- yt;
                cov_now = (eye(size(C,2)) - Kt1*C) * cov_next;
                mu_now  = mu_next + Kt1*rt;
                
                test(e,j) = yt;
                
                
                %2nd step update
                mu_vy           = mu_now(4); % mu_{t|t}^{V}
                sig_vy          = cov_now(4,4);
                
                myhat           = save_mu(1,index_ins);
                syhat           = save_cov(4,4,index_ins);
                
                mu_v2y          = mu_vy^2 + sig_vy; % mu_{t|t}^{V2}
                sig_v2y         = 2*sig_vy^2 + 4 * sig_vy * mu_vy^2;
                
                
                mu_v2           = myhat ;    % mu_{t|t-1}^{V2}
                cov_v2          = 3*syhat + 2*myhat^2;
                
                Kt2             = syhat/cov_v2;
                myhat           = myhat + Kt2*(mu_v2y - mu_v2); % mu_{t|t}^{V2hat}
                syhat           = syhat + (Kt2^2)*(sig_v2y - cov_v2);
                
                
                %sauvegarde
                
                %apply constraints on the deterioration speed
                if sqrt(myhat)<1
                    Dc=[1 ; 1];
                    d=[0.5 ; 1000000];
                    [myhat,syhat]=KFConstraintsHandlingSeq(myhat,syhat,Dc,d,1);
                end
                
                
                if (mu_now(2)+2*sqrt(cov_now(2,2)))>0 && j~=1
                    Dc=[0 1 0 0 ;0 1 0 0];
                    d=[-50;0];
                    [mu_now,cov_now]=KFConstraintsHandlingSeq(mu_now,cov_now,Dc,d,1);
                end
                
                
                
                save_mu(1,index_ins) = myhat;
                save_cov(:,:,index_ins) = cov_now(:,:);
                save_cov(4,4,index_ins) = syhat;
            end
            
            
        end
        save_epoch(i,:) = save_mu(1,:);
        save_cov(4,4,:) = save_mu(1,:);  %sig_v2
    end
    save_var_mu(e,:) = save_mu(1,:) ;
    save_var_sig(e,:) = save_cov(4,4,:);
    save_cov(4,4,:) = sd_sd_v0^2;
    
end

mu_v2hat = save_mu(1,:);

variance = mu_v2hat ;
sd = sqrt(variance);

end

