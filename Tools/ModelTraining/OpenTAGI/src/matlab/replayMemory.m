%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         replayMemory
% Description:  store and take samples for training RL agent
% Authors:      Luong-Ha Nguyen & James-A. Goulet
% Created:      September 16, 2020
% Updated:      October 27, 2020
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef replayMemory
    methods (Static)
        function [memory, count, filledIdx, fullMem] = push(memory, state, action, nextState, reward, finalState, count, maxMemory, fullMem, filledIdx)
            if count<= maxMemory                
                count = count+1;
                idx = count;
                if ~fullMem
                    filledIdx = count;
                end
            else
                filledIdx = maxMemory;
                count = 1;
                idx = count;
                fullMem = true;
            end           
            % Update memory
            memory(idx, 1) = {state};
            memory(idx, 2) = {action};
            memory(idx, 3) = {nextState};
            memory(idx, 4) = {reward};
            memory(idx, 5) = {finalState};
        end
        function [state, action, nextState, reward, finalState] = sample(memory, batchSize, filledIdx)
            idx        = randperm(filledIdx, batchSize); 
            samples    = memory(idx, :);
            state      = cat(1, samples{:, 1});
            action     = cat(1, samples{:, 2});
            nextState  = cat(1, samples{:, 3});
            reward     = cat(1, samples{:, 4});
            finalState = cat(1, samples{:, 5});
            finalState = finalState~=0;
        end
        function memory = init(maxMemory, dtype)
            memory = cell(maxMemory, 5);   
            memory(:, 1) = {zeros(1, 1, dtype)};
            memory(:, 2) = {zeros(1, 1, dtype)};
            memory(:, 3) = {zeros(1, 1, dtype)};
            memory(:, 4) = {zeros(1, 1, dtype)};
            memory(:, 5) = {zeros(1, 1, dtype)};
        end
        function [memory, count, fullMem, filledIdx] = push_V2(memory, state, action, nextState, nextAction, reward, finalState, count, maxMemory, fullMem, filledIdx)
            if count < maxMemory                
                count = count + 1;
                idx = count;
                if ~fullMem
                    filledIdx = count;
                end
            else
                filledIdx = maxMemory;              
                count     = 1;   
                idx       = count;
                fullMem   = true;
            end           
            % Update memory
            memory(idx, 1) = {state};
            memory(idx, 2) = {action};
            memory(idx, 3) = {nextState};
            memory(idx, 4) = {reward};
            memory(idx, 5) = {finalState};
            memory(idx, 6) = {nextAction};
        end
        function [state, action, nextState, nextAction, reward, finalState] = sample_V2(memory, batchSize, filledIdx)
            idx        = randperm(filledIdx, batchSize); 
            samples    = memory(idx, :);
            state      = cat(1, samples{:, 1});
            action     = cat(1, samples{:, 2});
            nextState  = cat(1, samples{:, 3});
            reward     = cat(1, samples{:, 4});
            finalState = cat(1, samples{:, 5});
            nextAction = cat(1, samples{:, 6});
        end
        function memory = init_V2(maxMemory, dtype)
            memory = cell(maxMemory, 5);   
            memory(:, 1) = {zeros(1, 1, dtype)};
            memory(:, 2) = {zeros(1, 1, dtype)};
            memory(:, 3) = {zeros(1, 1, dtype)};
            memory(:, 4) = {zeros(1, 1, dtype)};
            memory(:, 5) = {zeros(1, 1, dtype)};
            memory(:, 6) = {zeros(1, 1, dtype)};
        end
        function [memory, count, filledIdx] = pushAtari(memory, state, action, reward, finalState, count, maxMemory)
            if count<= maxMemory                
                count = count+1;
                idx = count;
                filledIdx = count;
            else
                filledIdx = maxMemory;
                count = 1;
                idx = count;
            end           
            % Update memory
            memory(idx, 1) = {state};
            memory(idx, 2) = {action};
            memory(idx, 3) = {reward};
            memory(idx, 4) = {finalState};
        end
        function [state, action, nextState, reward, finalState] = sampleAtari(memory, batchSize, numFrames, filledIdx)
            state      = cell(batchSize, 1);
            action     = zeros(batchSize, 1, 'like', memory{1, 2});
            nextState  = cell(batchSize, 1);
            reward     = zeros(batchSize, 1, 'like', memory{1, 3});
            finalState = zeros(batchSize, 1, 'like', memory{1, 4});
            loop   = 0;
            while loop < batchSize
                idx     = randperm(filledIdx-numFrames, 1);
                lastIdx = numFrames+1;
                samples = memory(idx:idx+lastIdx-1, :);
               finalStateloop = cat(1, samples{:, 4});              
                if all(finalStateloop(1:end-1)==0)
                    loop = loop + 1;
                    state{loop}      = cat(1, samples{1:lastIdx-1, 1});
                    action(loop)     = samples{lastIdx, 2};
                    nextState{loop}  = cat(1, samples{2:lastIdx, 1});                    
                    reward(loop)     = samples{lastIdx, 3};
                    finalState(loop) = samples{lastIdx, 4};
                end
            end
            state      = cat(1, state{:});
            nextState  = cat(1, nextState{:});
            finalState = finalState==1;
        end
        function memory = initAtari(maxMemory, dtype)
            memory = cell(maxMemory, 5);   
            memory(:, 1) = {zeros(1, 1, dtype)};
            memory(:, 2) = {zeros(1, 1, dtype)};
            memory(:, 3) = {zeros(1, 1, dtype)};
            memory(:, 4) = {zeros(1, 1, dtype)};
        end       
    end
end