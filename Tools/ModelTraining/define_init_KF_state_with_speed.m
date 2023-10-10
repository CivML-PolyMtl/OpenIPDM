function [init_x, init_V] = define_init_KF_state_with_speed(pred_vel, pred_vel_var, y, ave_init, Re, param)
    
%     max_obs = max(y,[],"all");
%     threshold = 0.95*max_obs;
%     mask = y(1,:,2) >= threshold;
%     pred_vel(mask) = 0;
%     pred_vel_var(mask) = mean(pred_vel_var(mask));
    init_x(2,:) = pred_vel; % output of ann_model
    init_V(2,2,:) = pred_vel_var; % output of ann_model
    
    init_x(1,:)=ave_init;
%     init_x(1,:)=gather(y(1,:,2));
%------
%     mask = y(1,:,3) > y(1,:,2);
%     init_x(1,mask) = (1-abs(mean(EngBiasData(:,2)))/75) .* max(y(1,mask,2:3),[], 3);
%     init_x(1,~mask) = y(1,~mask,2);
%     init_x(1,:) = max(y(1,:,2:3),[],3,'omitnan');
%------
%     init_x(1,:) = y(1,:,2) - EngBiasData(1,:,2);

%     init_x(1,:) = (1-abs(mean(EngBiasData(:,2)))/75) * max(gather(y(1,:,:)),[],3);
%     init_x(1,:)=(1 - abs(pred_vel.'/10)) .* max(gather(y_valid(1,:,:)),[],3);
%     init_x(1,:)=max((1 - abs(pred_vel.'/10)) .* max(gather(y_valid(1,:,:)),[],3), gather(y_valid(1,:,2)));
    init_x(3,:)=0;
    init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
    init_V(3,3,:)=param(5).^2;
end
