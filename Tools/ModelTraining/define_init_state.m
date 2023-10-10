classdef define_init_state
    methods(Static)
        function [init_x, init_V] = plain(y, Re, param)
            init_x(1,:)=gather(y(1,:,2));
            init_x(2,:)=0;
            init_x(3,:)=0;
        
            init_V=zeros(3,3,size(y,2),'gpuArray');
            init_V(3,3,:)=param(5).^2;
            init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
        end        
        function [init_x, init_V] = with_BNN(pred_vel, pred_vel_var, y, Re, param)
            init_x(1,:)=gather(y(1,:,2));
            init_x(2,:) = pred_vel;
            init_x(3,:)=0;
            
            init_V(1,1,:)=max(param(3).^2,Re(1,:,2));
            init_V(2,2,:) = pred_vel_var; 
            init_V(3,3,:)=param(5).^2;
        end
    end
end