%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         rltagi
% Description:  Selection actions, compute the Q values, normalization of
% states, output the video for Atari games.
% Authors:      Luong-Ha Nguyen & James-A. Goulet
% Created:      September 16, 2020
% Updated:      August 27, 2021
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef rltagi
    methods (Static) 
        % Action
        function action = getAction1net(net, theta, netStates, normStat,...
                maxIdx, state)
            if net.gpu
                state = gpuArray(state);
            end
            net.trainMode  = 0;
            netStates      = tagi.initializeInputs(netStates, state, [],...
                [], [], [], [], [], [], [], net.xsc);
            netStates      = tagi.feedForwardPass(net, theta, normStat,...
                netStates, maxIdx);
            [~, ~, mq, Sq] = tagi.extractStates(netStates);
            mql            = gather(mq{end});
            Sql            = gather(Sq{end});
            q              = normrnd(mql, sqrt(abs(Sql)));
            q(Sq{end}<0)   = -Inf;
            [~, action]    = max(q);
            if any(isnan(q))
                error('mean value of action distribution is nan')
            end
        end
        function action = getAction1netFullCov(net, theta, netStates, ...
                normStat, maxIdx, ms, Ssd, Ssf)
            if net.gpu
                ms = gpuArray(ms);
                Ssd = gpuArray(Ssd);
                Ssf = gpuArray(Ssf);
            end
            net.trainMode  = false;
            netStates      = tagi.initializeInputs(netStates, ms, Ssd, [],...
                [], [], [], [], [], [], net.xsc);
            netStates      = tagi.feedForwardPassFullCovCUDA(net, theta,...
                normStat, netStates, maxIdx, Ssf);
            [~, ~, mq, Sq] = tagi.extractStates(netStates);
            mql            = gather(mq{end});
            Sql            = gather(Sq{end});
            q              = normrnd(mql, sqrt(abs(Sql)));
            q(Sq{end}<0)   = -Inf;
            [~, action]    = max(q);
            if any(isnan(q))
                error('mean value of action distribution is nan')
            end
        end
        function action = getAction1netCUDA(net, theta, netStates, normStat,...
                maxIdx, state)
            if net.gpu
                state = gpuArray(state);
            end
            net.trainMode  = 0;
            netStates      = tagi.initializeInputs(netStates, state, [], [],...
                [], [], [], [], [], [], net.xsc);
            netStates      = tagi.feedForwardPassCUDA(net, theta, normStat,...
                netStates, maxIdx);
            [~, ~, mq, Sq] = tagi.extractStates(netStates);
            mql            = gather(mq{end});
            Sql            = gather(Sq{end});
            q              = normrnd(mql, sqrt(abs(Sql)));
            q(Sq{end}<0)   = -Inf;
            [~, action]    = max(q);
            if any(isnan(q))
                error('mean value of action distribution is nan')
            end
        end
        function [action, malA, SalA] = getContAction1net(netA, state)
            [thetaA, statesA, normStatA, maxIdxA] = network.extractNet(netA);         
            statesA = tagi.initializeInputs(statesA, state, [], [], [], [],...
                [], [], [], [], netA.xsc);
            statesA = tagi.feedForwardPass(netA, thetaA, normStatA, statesA,...
                maxIdxA);
            [~, ~, maA, SaA] = tagi.extractStates(statesA);
            malA = gather(maA{end});
            SalA = gather(SaA{end});
            action = normrnd(malA, sqrt(abs(SalA))); 
            action = cast(action, 'like', state);
        end
        
        % Rewards and Target value
        function [nMq, nSq, nAction] = nstepValue(netS, thetaS, statesS,...
                normStatS, maxIdxS, netT, thetaT, statesT, normStatT,...
                maxIdxT, state)
            statesS = tagi.initializeInputs(statesS, state, [], [], [], [],...
                [], [], [], [], netS.xsc); 
            statesS = tagi.feedForwardPass(netS, thetaS, normStatS, statesS,...
                maxIdxS);  
            [mzS, SzS, maS, SaS, JS, mdxsS, SdxsS,...
                mxsS, SxsS] = tagi.extractStates(statesS);
            
            statesT = tagi.initializeInputs(statesT, mzS{end}, SzS{end}, ...
                maS{end}, SaS{end}, JS{end}, mdxsS{end}, SdxsS{end}, ...
                mxsS{end}, SxsS{end}, netS.xsc);                        
            statesT = tagi.feedForwardPass(netT, thetaT, normStatT, statesT,...
                maxIdxT);  
               
            [~, ~, nextMq, nextSq] = tagi.extractStates(statesT);
            nextMq = reshape(nextMq{end}, ...
                [netT.ny, netT.batchSize * netT.repBatchSize]);
            nextSq = reshape(nextSq{end}, ...
                [netT.ny, netT.batchSize * netT.repBatchSize]);  
            [nMq, nSq, nAction] = rltagi.nextQvalues(nextMq, nextSq, netT.ny);
        end
        function [nMq, nSq, nAction] = nstepValue1net(netT, thetaT, statesT,...
                normStatT, maxIdxT, nextState)
            statesT = tagi.initializeInputs(statesT, nextState, [], [], [],...
                [], [], [], [], [], netT.xsc); 
            statesT = tagi.feedForwardPassCUDA(netT, thetaT, normStatT, ...
                statesT, maxIdxT);                
            [~, ~, nextMq, nextSq] = tagi.extractStates(statesT);
            nextMq = reshape(gather(nextMq{end}), ...
                [netT.ny, netT.batchSize * netT.repBatchSize]);
            nextSq = reshape(gather(nextSq{end}), ...
                [netT.ny, netT.batchSize * netT.repBatchSize]);  
            [nMq, nSq, nAction] = rltagi.nextQvalues(nextMq, nextSq, netT.ny);
        end
        function [nMq, nSq] = nstepValue_V2(netQ, netA, nextState,...
                nextAction)
            actFunIdx = 1;
            [thetaQ, statesQ, normStatQ, maxIdxQ] = network.extractNet(netQ);          
            [malA] = act.meanVar(nextAction, nextAction, ...
                zeros(size(nextAction)), actFunIdx, netQ.actBound(1), ...
                netQ.batchSize, netQ.repBatchSize, netQ.gpu) ;           
            ma0Q = tagi.attachMeanVar(nextState, malA, netQ.nx - netA.ny, ...
                netA.ny, netQ.batchSize, netQ.repBatchSize);
            
            statesQ = tagi.initializeInputs(statesQ, ma0Q, [], [], [], [], ...
                [], [], [], [], netQ.xsc);
            statesQ = tagi.feedForwardPassCUDA(netQ, thetaQ, normStatQ, ...
                statesQ, maxIdxQ);
            [~, ~, nextMq, nextSq] = tagi.extractStates(statesQ);
            nMq = nextMq{end};
            nSq = nextSq{end};            
        end
        function [malQ, SalQ] = nstepValueBatch(netQ, netA, state, action)
            [thetaQ, statesQ, normStatQ, maxIdxQ] = network.extractNet(netQ);                                           
            % Get Q values          
            [malA] = act.meanVar(action, action, zeros(size(action)), ...
                netQ.actFunIdx(1), netQ.actBound(1), netQ.batchSize, ...
                netQ.repBatchSize, netQ.gpu) ;  
            ma0Q = tagi.attachMeanVar(state, malA, netQ.nx - netA.ny, ...
                netA.ny, netQ.batchSize, netQ.repBatchSize);
            
            statesQ  = tagi.initializeInputs(statesQ, ma0Q, [], [], [], [],...
                [], [], [], [], netQ.xsc);
            statesQ = tagi.feedForwardPass(netQ, thetaQ, normStatQ, statesQ,...
                maxIdxQ);
            [~, ~, maQ, SaQ] = tagi.extractStates(statesQ);
            malQ = maQ{end};
            SalQ = SaQ{end};                     
        end
        function [nMq, nSq] = nstepValue_V3(netQ, netA, nextState, ...
                nextAction)
            onePadS = ones(size(nextState));
            [thetaQ, statesQ, normStatQ, maxIdxQ] = network.extractNet(netQ);          
            mz0Q = tagi.attachMeanVar(nextState, nextAction, ...
                netQ.nx - netA.ny, netA.ny, netQ.batchSize, netQ.repBatchSize);
            [malA,~,JlA] = act.meanVar(nextAction, nextAction, ...
                zeros(size(nextAction)), netQ.actFunIdx(1), netQ.actBound(1),...
                netQ.batchSize, netQ.repBatchSize, netQ.gpu) ;  
            
            ma0Q = tagi.attachMeanVar(nextState, malA, netQ.nx - netA.ny, ...
                netA.ny, netQ.batchSize, netQ.repBatchSize);
            J0Q = tagi.attachMeanVar(onePadS, JlA, netQ.nx-netA.ny, netA.ny,...
                netQ.batchSize, netQ.repBatchSize);
            
            statesQ = tagi.initializeInputs(statesQ, mz0Q, [], ma0Q, [],...
                J0Q, [], [], [], [], netQ.xsc);
            [statesQ] = tagi.feedForwardPass(netQ, thetaQ, normStatQ, ...
                statesQ, maxIdxQ);
            [~, ~, nextMq, nextSq] = tagi.extractStates(statesQ);
            nMq = nextMq{end-1};
            nSq = nextSq{end-1};            
        end
        function [nMq, nSq] = nstepLastStep(netQ, netS, nextState, nextAction)
            zeroPadA = ones(size(nextAction));
            [thetaS, statesS, normStatS, maxIdxS] = network.extractNet(netS);
            [thetaQ, statesQ, normStatQ, maxIdxQ] = network.extractNet(netQ);
            
            % Get Q values
            statesS   = tagi.initializeInputs(statesS, nextState, [], [], ...
                [], [], [], [], [], [], netS.xsc);
            [statesS] = tagi.feedForwardPass(netS, thetaS, normStatS, ...
                statesS, maxIdxS);            
            
            [mzS, SzS, maS, SaS, JS] = tagi.extractStates(statesS);
            [malA, ~, JlA] = act.meanVar(nextAction, nextAction, zeroPadA,...
                netQ.actFunIdx(1), netQ.actBound(1), netQ.batchSize, ...
                netQ.repBatchSize, netQ.gpu) ; 
            
            mz0Q = tagi.attachMeanVar(mzS{end}, nextAction, netS.ny, ...
                netQ.nx - netS.ny, netQ.batchSize, netQ.repBatchSize);
            Sz0Q = tagi.attachMeanVar(SzS{end}, zeroPadA, netS.ny, ...
                netQ.nx - netS.ny, netQ.batchSize, netQ.repBatchSize);
            ma0Q = tagi.attachMeanVar(maS{end}, malA, netS.ny, ...
                netQ.nx - netS.ny, netQ.batchSize, netQ.repBatchSize);
            Sa0Q = tagi.attachMeanVar(SaS{end}, zeroPadA, netS.ny, ...
                netQ.nx - netS.ny, netQ.batchSize, netQ.repBatchSize);
            J0Q  = tagi.attachMeanVar(JS{end}, JlA, netS.ny, ...
                netQ.nx - netS.ny, netQ.batchSize, netQ.repBatchSize);
            
            statesQ   = tagi.initializeInputs(statesQ, mz0Q, Sz0Q, ma0Q, ...
                Sa0Q, J0Q, [], [], [], [], netQ.xsc);
            [statesQ] = tagi.feedForwardPass(netQ, thetaQ, normStatQ, ...
                statesQ, maxIdxQ);
            [~, ~, nextMq, nextSq] = tagi.extractStates(statesQ);
            nMq = nextMq{end};
            nSq = nextSq{end};            
        end

        function [Sy] = nstepSv(Sv, dones, gamma, nsteps)
            Sy = zeros(nsteps, 1, 'like', Sv); 
            SR = Sv;
            for n = nsteps:-1:1
                if dones(n + 1)
                    SR = Sv;
                end
                SR    = (gamma.^2) * SR;
                Sy(n) = SR;
            end           
        end
        function [Sy] = nstepSv_V2(Sv, dones, gamma, nsteps)
            Sy = zeros(nsteps, 1, 'like', Sv); 
            SR = Sv;
            for n = 1:nsteps
                if dones(n + 1)
                    SR = Sv;
                end
                SR    = SR ./ (gamma.^2) ;
                Sy(n) = SR;
            end           
        end  
        function [Sy] = nstepSv_V3(Sv, dones, gamma, nsteps)
            Sy = zeros(nsteps, 1, 'like', Sv); 
            countDone = 0;
            loop = 0;
            for n = nsteps:-1:1                
                if dones(n + 1)
                    countDone = countDone + 1;
                end
                if countDone >= 1
                    break;
                end
                Sy(n) = 1;
                loop  = loop + 1;
            end 
            Syloop = zeros(loop, 1, 'like', Sv); 
            for i = 1:loop
                if i>1                   
                    Sv = Sv ./ (gamma.^2);
                end
                Syloop(i) = Sv;
            end
            Sy(Sy==1) = Syloop;
        end        
        
        function [nextMq, nextSq, idxMax] = nextQvalues(mq, Sq, ny)
            nextq       = normrnd(mq, sqrt(abs(Sq)));
            nextq(Sq<0) = -Inf;
            [~, idxMax] = max(nextq);
            idxMax      = idxMax + colon(0, ny, numel(mq)-ny);
            nextMq      = mq(idxMax)';
            nextSq      = Sq(idxMax)';
            idxMax      = idxMax(:);
        end
        function discountRew = discountReward(rewards, gamma, totalSteps)
            discountRew = zeros(totalSteps, 1, 'like', rewards);
            R = 0;
            for n = totalSteps:-1:1
                R = gamma*R + rewards(n);
                discountRew(n) = R;
            end         
        end
        function [nextdMq, nextdSq, CnextdQy] = discountValue(nextMq, nextSq,...
                gamma, totalSteps)
            n = totalSteps:-1:1;
            n = n(:);
            nextdMq  = (gamma.^(n)).*nextMq;
            nextdSq  = (gamma.^(2*n)).*(nextSq);
            CnextdQy = -(gamma.^(n)).*nextSq;
        end
        function [nextdMq, nextdSq, CnextdQy] = discountValueCont(nextMq,...
                nextSq, gamma, totalSteps)
            N = size(nextMq, 1);
            nextdMq  = zeros(N, totalSteps, 'like', nextMq);
            nextdSq  = zeros(N, totalSteps, 'like', nextSq);
            CnextdQy = zeros(N, totalSteps, 'like', nextSq);
            for n = totalSteps:-1:1
                nextdMq(:, n)  = (gamma^(totalSteps-n+1))*nextMq;
                nextdSq(:, n)  = (gamma^(2*(totalSteps-n+1)))*nextSq;
                CnextdQy(:, n) = -(gamma^(totalSteps-n+1))*nextSq;
            end
        end
        function x = clamp(x, xmin, xmax) 
            if any(x<xmin)
                x(x<xmin) = xmin(x<xmin);
            end
            if any(x>xmax)
                x(x>xmax) = xmax(x>xmax);
            end
        end        
        
        % Normalization state
        function [newMs, newSs, newCount] = runningMeanStd(state, ms, Ss, ...
                count, batchCount, numStates)
            % This code is based on openai's code
            batchCount= cast(batchCount, 'like', ms);
            state     = cat(1, state{:, 1});
            state     = reshape(state, [numStates, batchCount]);
            batchMean = mean(state, 2);
            batchStd  = std(state, 0, 2);
            batchStd(batchStd==0) = 1;
            batchVar  = batchStd.^2;
            
            delta    = batchMean - ms;
            totCount = count + batchCount;
            newMs    = ms + delta * (batchCount / totCount);
            M2a      = Ss * count;
            M2b      = batchVar * batchCount;
            M2       = M2a + M2b + (delta.^2) * count * (batchCount / totCount);
            newSs    = M2 / totCount;
            newCount = totCount;
        end
        function [newMs, newSs, newCount] = runningMeanStd_V2(x, mx, Sx, ...
                count, batchCount)
            % This code is based on openai's code
            batchCount = cast(batchCount, 'like', mx);
            batchMean  = mean(x, 2);
            if batchCount==1
                batchStd = 0;
            else
                batchStd = std(x, 0, 2);
            end
            batchVar   = batchStd.^2;
            
            delta    = batchMean - mx;
            totCount = count + batchCount;
            newMs    = mx + delta * (batchCount / totCount);
            M2a      = Sx * count;
            M2b      = batchVar * batchCount;
            M2       = M2a + M2b + (delta.^2) * count * (batchCount / totCount);
            newSs    = M2 / totCount;
            newCount = totCount;
        end
        function [newMs, newSs, newCount] = runningMeanStd_V3(x, mx, Sx, ...
                count, batchCount)
            % This code is based on openai's code
            batchCount = cast(batchCount, 'like', mx);
            batchMean  = mean(x, 1);
            if batchCount==1
                batchStd = 0;
            else
                batchStd = std(x,1);
            end
            batchVar   = batchStd.^2;
            
            delta    = batchMean - mx;
            totCount = count + batchCount;
            newMs    = mx + delta * (batchCount / totCount);
            M2a      = Sx * count;
            M2b      = batchVar * batchCount;
            M2       = M2a + M2b + (delta.^2) * count * (batchCount / totCount);
%             newSs    = M2 / (totCount);
            newSs    = M2 / (totCount - 1);
            newCount = totCount;
        end
        function [newMs, newSs, newCount] = runningMeanStdDist(x, Sx, rmx, ...
                rSx, count, batchCount)
            % This code is based on openai's code
            batchCount = cast(batchCount, 'like', rmx);
            batchMean  = mean(x, 1);
            batchVar   = 1 / batchCount * (sum(Sx, 1)+ sum((x - batchMean).^2, 1));
            
            delta    = batchMean - rmx;
            totCount = count + batchCount;
            newMs    = rmx + delta * (batchCount / totCount);
            M2a      = rSx * count;
            M2b      = batchVar * batchCount;
            M2       = M2a + M2b + (delta.^2) * count * (batchCount / totCount);
            newSs    = M2 / totCount;
            newCount = totCount;
        end       
        
        % Video
        function video(images, vname)
            writerObj           = VideoWriter(vname);
            writerObj.FrameRate = 100;
            open(writerObj);
            for i = 1:length(images)
                frame = im2frame(uint8(images{i}));
                writeVideo(writerObj, frame);
            end
            close(writerObj);
        end
        function video_V2(frames)
           frames = cat(4, frames{:});
           frames = frames(:,:,1,:);
           mov    = immovie(uint8(frames), colormap);
           implay(mov)
        end
    end
end