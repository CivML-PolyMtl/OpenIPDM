%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         wrappers
% Description:  preprocess env for atari games 
% Authors:      Luong-Ha Nguyen & James-A. Goulet
% Created:      September 16, 2020
% Updated:      August 19, 2021
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef wrappers
    methods (Static)
        function resizedImg = resizeImage(image, method)
            if strcmp(method, 'scale')
                resizedImg = imresize(image, [84, 84], 'bilinear', 'Antialiasing', false);
                resizedImg = single(resizedImg);
            elseif strcmp(method, 'crop')
                resizedImg = imresize(image, [160, 84], 'bilinear', 'Antialiasing', false);
                resizedImg = resizedImg(18:101, :); 
            else
                resizedImg = image;
            end
        end
        function img = atariPreprocessing(img, method)
            img = rgb2gray(img);
            img = wrappers.resizeImage(img, method);
        end
        function [state, obs, totalReward, done, info] = maxAndSkipEnv(env, lives, action, skip)
            obs          = cell(skip, 1);
            obsGrayscale = cell(skip, 1);
            totalReward  = 0;
            for s = 1:skip
                [img, reward, done, info] = env.step(action);
                newlives = info;
                totalReward = totalReward + reward;
                obs{s} = img;
                
                obsGrayscale{s} = img;
                if done||(newlives<lives&&newlives>0)
                    obs = obs(1:s);
                    obsGrayscale = obsGrayscale(1:s);
                    break
                end               
            end
            if size(obsGrayscale, 1)~=1
                state = max(obsGrayscale{end}, obsGrayscale{end-1});
            else 
                state = obsGrayscale{1};
            end
        end
        function [state, rewardclip, reward, obs, done, info] = stackFrames(env, action, state, lives, numFrames, skip, resizeMethod)
            % Initialize state
            state(:,:,1:numFrames-1) = state(:,:,2:numFrames);
            
            % Skip frames
            [nextState, obs, reward, done, info] = wrappers.maxAndSkipEnv(env, lives, action, skip);
            nextState  = rgb2gray(nextState);
            resizedImg = wrappers.resizeImage(nextState, resizeMethod);
            state(:, :, numFrames) = resizedImg ./ 255;

            % Clip reward 
            rewardclip = sign(reward);
        end
        function [state, rewardclip, reward, obs, done, info] = getFrames(env, action, lives, skip, resizeMethod)
            [nextState, obs, reward, done, info] = wrappers.maxAndSkipEnv(env, lives, action, skip);
            nextState  = rgb2gray(nextState);
            resizedImg = wrappers.resizeImage(nextState, resizeMethod);
            state = resizedImg ./ 255;
            
            % Clip reward 
            rewardclip = sign(reward);
        end
        function [obs, reward, done, info] = noopResetEnv(env)
            obs = env.reset();
            noops = randi(30, 1);            
            for i = 1:noops
                [obs, reward, done, info] = env.step(int16(0));
                if done
                    obs = env.reset();
                    error('there is a problem in noop start')                  
                end
            end            
        end
        function obs = fireResetEnv(env)
            obs = env.reset();
            [obs, ~, done] = env.step(int16(1));
            if done
                obs = env.reset();
            end
            [obs, ~, done] = env.step(int16(2));
            if done
                obs = env.reset();
            end
        end
        function obs = fireNoopResetEnv(env)
            obs = wrappers.fireResetEnv(env);
            noops = randi(30, 1);            
            for i = 1:noops
                [obs, ~, done] = env.step(int16(0));
                if done == 1
                    obs = wrappers.fireResetEnv(env);
                end
            end
        end
        function [obs, action] = noopFireResetEnv(env)      
            run = true;
            while run
                [~, done] = wrappers.noopResetEnv(env);
                if done
                    error('there is a problem in noopfire start')
                end
                [obs, ~, done] = env.step(1);
                action = 2;
                if~done
                    [obs, ~, done] = env.step(2);
                    action = 3;
                end
                if~done; break; end
            end
        end
        function [state, obs, info] = noopFireResetEnv_V3(env, skip, resizeMethod) 
            % Initialization
            obs = cell(skip*2, 1);
            
            % Run loop
            run = true;
            while run                
                % Noop reset
                [obs, ~, done, info] = wrappers.noopResetEnv(env);
                lives = info;                
                if done
                    error('there is a problem in noopfire start')
                end
                
                % Fire reset
                action = 1;
                [state, obs, ~, done, info] = wrappers.maxAndSkipEnv(env, lives, action, skip);
                
                frameCount1 = size(obs, 1);
                obs(1:frameCount1, 1) = obs;
                if~done; break; end                
            end
            obs   = obs(1:frameCount1);
            state = rgb2gray(state);
            state = wrappers.resizeImage(state, resizeMethod);
            state = state./255;
        end
        function [state, obs, reward, done, info, env] = episodicLifeResetEnv(env, lives, skip, resizeMethod)            
            % Initialization
            obs = cell(skip*3, 1);
            
            % Noop actions
            action = 0;   
            [~, obs0, ~, done, info] = wrappers.maxAndSkipEnv(env, lives, action, skip);
            lives = info;            
            frameCount0 = size(obs0, 1);
            obs(1:frameCount0, 1) = obs0;
            if done
                error('There is a problem with episodic life');
            end
            
            % Fire reset
            action = 1;
            [state, obs1, reward, done, info] = wrappers.maxAndSkipEnv(env, lives, action, skip);
                     
            frameCount1 = size(obs1, 1) + frameCount0;
            obs(frameCount0+1:frameCount1, 1) = obs1;
            obs   = obs(1:frameCount1);
            state = rgb2gray(state);
            state = wrappers.resizeImage(state, resizeMethod);
            state = state./255;
        end
    end
end