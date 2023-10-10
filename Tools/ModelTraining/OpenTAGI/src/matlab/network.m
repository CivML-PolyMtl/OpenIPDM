%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         network
% Description:  Build networks relating to each task (task.m)
% Authors:      Luong-Ha Nguyen & James-A. Goulet & Ali Fakhri 
% Created:      July 02, 2020
% Updated:      September 15, 2023 | Removed functions not used by OpenIPDM
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca & afakhri@pm.me
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef network
    methods (Static)
        % Regression 
        function [theta, normStat, zl, Szl, sv] = regression(net, theta,...
                normStat, states, maxIdx, x, y, Sy)
            % Initialization
            numObs = size(x, 1);
            numDataPerBatch = net.repBatchSize * net.batchSize;
            zl  = zeros(numObs, net.ny, net.dtype);
            Szl = zeros(numObs, net.ny, net.dtype);
            % Loop
            loop = 0;
            for i = 1 : numDataPerBatch : numObs
                loop     = loop + 1; 
                if numDataPerBatch==1
                    idxBatch = i : i + net.batchSize - 1;
                else
                    if numObs - i >= numDataPerBatch
                        idxBatch = i : i + numDataPerBatch - 1;
                    else
                        idxBatch = [i : numObs, randperm(i - 1, ...
                            numDataPerBatch - numObs + i - 1)];
                    end
                end
                % Covariate
                xloop  = reshape(x(idxBatch, :)', ...
                    [net.batchSize * net.nx, net.repBatchSize]);
                states = tagi.initializeInputs(states, xloop, [], [], [],...
                    [], [], [], [], [], net.xsc);
                % Training
                if net.trainMode                 
                    % Observation
                    yloop = reshape(y(idxBatch, :)', ...
                        [net.batchSize * net.nl, net.repBatchSize]);
                    if net.imperfect_obs
                        Syloop = reshape(Sy(idxBatch, :)', ...
                            [net.batchSize * net.nl, net.repBatchSize]);
                    else
                        Syloop = [];
                    end
                    [states, normStat, maxIdx] = tagi.feedForwardPass(net,...
                        theta, normStat, states, maxIdx);
                    [deltaM, deltaS,deltaMx, deltaSx,...
                        ~, ~, sv] = tagi.hiddenStateBackwardPass(net, theta,...
                        normStat, states, yloop, Syloop, [], maxIdx);
                    net.sv = sv;
                    dtheta = tagi.parameterBackwardPass(net, theta, normStat,...
                        states, deltaM, deltaS, deltaMx, deltaSx);
                    % ------------------------------------constraint the update to parameters
                    % get min between abs of dtheta or 0.5 * theta
                    get_min = @(x, y) sign(x) .* min(abs(x), 0.5.*abs(y));
                    dtheta(1:4,:) = cellfun(get_min, dtheta(1:4,:), theta(1:4,:), 'UniformOutput', false);

                    % constraint to [-1,1] in case theta has large values..
                    max_constraint = @(x) min(x, 1);
                    dtheta(1:4,:) = cellfun(max_constraint, dtheta(1:4,:), 'UniformOutput', false);
                    min_constraint = @(x) max(x,-1);
                    dtheta(1:4,:) = cellfun(min_constraint, dtheta(1:4,:), 'UniformOutput', false);
                    % ------------------------------------constraint the update to parameters
                    
                    theta  = tagi.globalParameterUpdate(theta, dtheta,...
                        net.gpu);  
                    [~, ~, ma, Sa]   = tagi.extractStates(states);
                % Testing    
                else 
                    [states, normStat, maxIdx] = tagi.feedForwardPass(net,...
                        theta, normStat, states, maxIdx);
                    [~, ~, ma, Sa]   = tagi.extractStates(states);
                    zl(idxBatch, :)  = gather(reshape(ma{end}, ...
                        [net.ny, numDataPerBatch])');
                    Szl(idxBatch, :) = gather(reshape(Sa{end}, ...
                        [net.ny, numDataPerBatch])');
                    sv = net.sv;
                end 
                zl(idxBatch, :)  = gather(reshape(ma{end}, ...
                    [net.ny, numDataPerBatch])');
                Szl(idxBatch, :) = gather(reshape(Sa{end}, ...
                    [net.ny, numDataPerBatch])');
            end
        end
               
        % Initialization
        function [net, states, maxIdx, netInfo] = initialization(net)
            % Build indices
            net = indices.initialization(net);
            net = indices.layerEncoder(net);
            net = indices.parameters(net);
            if net.cuda 
                net = indices.covarianceCUDA(net);
            else
                net = indices.covariance(net);
            end
            netInfo = indices.savedInfo(net);
            % States
            states = tagi.initializeStates(net.nodes, net.batchSize, ...
                net.repBatchSize, net.xsc, net.dtype, net.gpu);
            maxIdx = tagi.initializeMaxPoolingIndices(net.nodes, net.layer,...
                net.layerEncoder, net.batchSize, net.repBatchSize, ...
                net.dtype, net.gpu);                        
        end
        function [theta, states, normStat, maxIdx] = extractNet(net)
            theta    = net.theta;
            states   = net.states;
            normStat = net.normStat;
            maxIdx   = net.maxIdx;
        end
        function net = compressNet(net, theta, states, normStat, maxIdx)
            net.theta    = theta;
            net.states   = states;
            net.normStat = normStat;
            net.maxIdx   = maxIdx;
        end
    end
end