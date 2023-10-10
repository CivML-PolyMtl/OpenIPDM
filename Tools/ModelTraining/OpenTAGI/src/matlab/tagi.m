%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         tagi
% Authors:      Luong-Ha Nguyen & James-A. Goulet & Ali Fakhri 
% Created:      November 03, 2019
% Updated:      September 15, 2023 | Removed functions not used by OpenIPDM
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca & afakhri@pm.me
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Description: Tractable Approximate Gaussian Inference (TAGI). 
% This function contains the core of TAGI that includes
% - forward uncertain propagation
% - backward update i.e., inference.
% Further details for each steps, please check out this following paper
% Goulet, Nguyen, and Amiri, (2020): Tractable Approximate Gaussian
% Inference for Bayesian Neural Networks. 
% These two steps are implemented for different types of neural netwworks
% such as FNNs, CNNs, pooling, normalization etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 classdef tagi
    methods(Static) 
        % Feedforward
        function [states, normStat, maxIdx, mda, Sda] = feedForwardPass(net,...
                theta, normStat, states, maxIdx)
            % Initialization
            [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = tagi.extractParameters(theta);
            [mz, Sz, ma, Sa, J, mdxs, Sdxs,...
                mxs, Sxs] = tagi.extractStates(states);
            [mra, Sra] = tagi.extractNormStat(normStat);
            numLayers  = length(net.nodes);
            actFunIdx  = net.actFunIdx;
            actBound   = net.actBound;
            layer      = net.layer;
            imgW       = cast(net.imgW, net.dtype);
            imgH       = cast(net.imgH, net.dtype);
            filter     = cast(net.filter, net.dtype);
            kernelSize = cast(net.kernelSize, net.dtype);
            B          = cast(net.batchSize, net.dtype);
            rB         = cast(net.repBatchSize, net.dtype);
            nodes      = cast(net.nodes, net.dtype);
            epsilon    = net.epsilon;
            mhat       = cell(numLayers, 1);
            Shat       = cell(numLayers, 1);
            numParamsPerlayer_2 = net.numParamsPerlayer_2;
            smlIdx     = net.similarIdx;
            
            % Derivative
            mda    = cell(numLayers, 1);
            Sda    = cell(numLayers, 1);
            mda{1} = ones(size(mz{1}), 'like', mz{1});
            Sda{1} = zeros(size(Sz{1}), 'like', Sz{1});
            
            % Hidden Layers
            for j = 2:numLayers
                idxw = (numParamsPerlayer_2(1, j-1)+1):numParamsPerlayer_2(1, j);
                idxb = (numParamsPerlayer_2(2, j-1)+1):numParamsPerlayer_2(2, j);
                % Max pooling
                if layer(j) == net.layerEncoder.mp
                    maPool = normrnd(gather(ma{j-1}), sqrt(abs(gather(Sa{j-1}))));
                    if net.padding(j-1) ~= 0
                        maPool = vertcat(maPool,...
                            -Inf * ones(1, size(maPool, 2), 'like', maPool));
                    end
                    maPool(Sa{j-1}<=0) = -Inf;
                    [mz{j}, Sz{j}, maxIdx{j}] = tagi.mpMeanVar(mz{j}, Sz{j},...
                        maPool, ma{j-1}, Sa{j-1}, net.idxPooling{j-1},...
                        maxIdx{j}, rB, net.gpu);
                    
                    % Average pooling
                elseif layer(j) == net.layerEncoder.ap
                    [mz{j}, Sz{j}] = tagi.apMeanVar(mz{j}, Sz{j}, ma{j-1},...
                        Sa{j-1}, net.idxPooling{j-1}, net.padding(j-1), rB);
                    
                    % Normalization
                elseif layer(j) == net.layerEncoder.ln...
                        || layer(j) == net.layerEncoder.bn
                    
                    if net.trainMode
                        [mhat{j-1}, Shat{j-1}] = tagi.pMeanVar(ma{j-1},...
                            Sa{j-1}, nodes(j-1), imgW(j-1), imgH(j-1),...
                            filter(j-1), B, rB, layer(j-1), layer(j),...
                            net.layerEncoder);
                        % Running average for mean and variance
                        mra{j-1} = net.normMomentum*mra{j-1} +...
                            (1 - net.normMomentum) * mhat{j-1};
                        Sra{j-1} = net.normMomentum*Sra{j-1} +...
                            (1 - net.normMomentum) * Shat{j-1};
                    end
                    mhatD = tagi.distributeNormMeanVar(mra{j-1}, nodes(j-1),...
                        imgW(j-1), imgH(j-1), filter(j-1), B, rB,...
                        layer(j-1), layer(j), net.layerEncoder);
                    ShatD = tagi.distributeNormMeanVar(Sra{j-1}, nodes(j-1),...
                        imgW(j-1), imgH(j-1), filter(j-1), B, rB,...
                        layer(j-1), layer(j), net.layerEncoder);
                    if layer(j-1) == net.layerEncoder.fc
                        [mz{j}, Sz{j}] = tagi.fcNormMeanVar(mz{j}, Sz{j},...
                            mw(idxw), Sw(idxw), mb(idxb), Sb(idxb), ma{j-1},...
                            Sa{j-1}, mhatD, ShatD, epsilon, B, rB, net.gpu);
                    elseif layer(j-1) == net.layerEncoder.conv...
                            || layer(j-1) == net.layerEncoder.tconv
                        [mz{j}, Sz{j}] = tagi.convNormMeanVar(mz{j}, Sz{j},...
                            mw(idxw), Sw(idxw), mb(idxb), Sb(idxb), ma{j-1},...
                            Sa{j-1}, mhatD, ShatD, epsilon, imgH(j-1),...
                            imgH(j-1), filter(j-1), B, rB, net.gpu);
                    end
                    
                    % Convolutional
                elseif layer(j) == net.layerEncoder.conv
                    if B==1&&rB==1
                        [mz{j}, Sz{j}] = tagi.convMeanVarB1(mw(idxw),...
                            Sw(idxw), mb(idxb), Sb(idxb), ma{j-1}, Sa{j-1},...
                            net.idxFmwa(smlIdx(j-1), :), kernelSize(j-1),...
                            filter(j-1), imgW(j), imgH(j), filter(j),...
                            net.padding(j-1), net.gpu);
                    else
                        [mz{j}, Sz{j}] = tagi.convMeanVar(mz{j}, Sz{j},...
                            mw(idxw), Sw(idxw), mb(idxb), Sb(idxb), ma{j-1},...
                            Sa{j-1}, net.idxFmwa(smlIdx(j-1), :),...
                            kernelSize(j-1), filter(j-1), imgW(j), imgH(j),...
                            filter(j), B, rB, net.padding(j-1), net.gpu);
                    end
                    % Transposed convolutional
                elseif layer(j) == net.layerEncoder.tconv
                    [mz{j}, Sz{j}] = tagi.tconvMeanVar(mz{j}, Sz{j},...
                        mw(idxw), Sw(idxw), mb(idxb), Sb(idxb), ma{j-1},...
                        Sa{j-1}, net.idxFmwa(j-1, :), imgW(j), imgH(j),...
                        filter(j), B, rB, net.gpu);
                    
                    % Full-connected
                elseif layer(j) == net.layerEncoder.fc
                    [mz{j}, Sz{j}] = tagi.fcMeanVar(mz{j}, Sz{j}, mw(idxw),...
                        Sw(idxw), mb(idxb), Sb(idxb), ma{j-1}, Sa{j-1},...
                        nodes(j-1), nodes(j), B, rB, net.gpu);
                end
                
                % Shortcut connection for residual networks
                if net.xsc(j)~=0 && (net.filter(net.xsc(j))~=net.filter(j)...
                        ||net.imgW(net.xsc(j))~=net.imgW(j))
                    idxXsc = net.xsc(j);
                    idxwx = (numParamsPerlayer_2(3, idxXsc)+1):numParamsPerlayer_2(3, idxXsc+1);
                    idxbx = (numParamsPerlayer_2(4, idxXsc)+1):numParamsPerlayer_2(4, idxXsc+1);
                    [mxs{j}, Sxs{j}] = tagi.convMeanVar(mxs{j}, Sxs{j},...
                        mwx(idxwx), Swx(idxwx), mbx(idxbx), Sbx(idxbx),...
                        ma{idxXsc}, Sa{idxXsc}, net.idxFmwaXsc(idxXsc, :),...
                        1, filter(idxXsc), imgW(j), imgH(j), filter(j), B,...
                        rB, net.paddingXsc(idxXsc), net.gpu);
                    % Save convolutional hidden state before adding x
                    % shortcut
                    mdxs{j} = mz{j};
                    Sdxs{j} = Sz{j};
                    [mz{j}, Sz{j}] = arrayfun(@twoPlus, mz{j}, Sz{j},...
                        mxs{j}, Sxs{j});
                elseif net.xsc(j)~=0&&(net.filter(net.xsc(j))==net.filter(j)...
                        ||net.imgW(net.xsc(j))~=net.imgW(j))
                    mxs{j}  = mz{net.xsc(j)};
                    Sxs{j}  = Sz{net.xsc(j)};
                    mdxs{j} = mz{j};
                    Sdxs{j} = Sz{j};
                    [mz{j}, Sz{j}] = arrayfun(@twoPlus, mz{j}, Sz{j},...
                        mxs{j}, Sxs{j});
                end
                
                % Activation
                if actFunIdx(j)~=0
                    if net.split_v2_act && j == numLayers
                        [ma{j}(1:net.nl), Sa{j}(1:net.nl), J{j}(1:net.nl)] = act.meanVar(mz{j}(1:net.nl), mz{j}(1:net.nl), Sz{j}(1:net.nl),...
                        actFunIdx(j), actBound(j), B, rB, net.gpu);
                        ma{j}(net.nl+1:end) = mz{j}(net.nl+1:end);
                        Sa{j}(net.nl+1:end) = Sz{j}(net.nl+1:end);
                        J{j}(net.nl+1:end)  = ones(size(mz{j}(net.nl+1:end)), 'like', mz{j}(net.nl+1:end));
                    else    
                        [ma{j}, Sa{j}, J{j}] = act.meanVar(mz{j}, mz{j}, Sz{j},...
                        actFunIdx(j), actBound(j), B, rB, net.gpu);
                    end
                else
                    ma{j} = mz{j};
                    Sa{j} = Sz{j};
                    J{j}  = ones(size(mz{j}), 'like', mz{j});
                end
                
                % Derivative for FC
                if net.collectDev&&actFunIdx(j)~=0
                    [mda{j}, Sda{j}] = act.meanVarDev(mz{j}, Sz{j},...
                        actFunIdx(j), actBound(j));
                end
            end

            [mz, Sz, ma, Sa] = cap_mean_var(mz, Sz, ma, Sa);
            normStat = tagi.compressNormStat(mra, Sra);
            states   = tagi.compressStates(mz, Sz, ma, Sa, J, mdxs, Sdxs,...
                mxs, Sxs);
        end
        
        % Inference 
        function [deltaM, deltaS, deltaMx, deltaSx,...
                deltaMz0, deltaSz0, sv] = hiddenStateBackwardPass(net,...
                theta, normStat, states, y, Sy, udIdx, maxIdx)
            % Initialization
            [mw, ~, ~, ~, mwx] = tagi.extractParameters(theta);
            [mz, Sz, ma, Sa, J, mdxs, Sdxs,...
                ~, Sxs] = tagi.extractStates(states);
            [~, Sra] = tagi.extractNormStat(normStat);
            numLayers  = length(net.nodes);
            imgW       = cast(net.imgW, net.dtype);
            imgH       = cast(net.imgH, net.dtype);
            filter     = cast(net.filter, net.dtype);
            kernelSize = cast(net.kernelSize, net.dtype);
            stride     = cast(net.stride, net.dtype);
            B          = cast(net.batchSize, net.dtype);
            rB         = cast(net.repBatchSize, net.dtype);
            nodes      = cast(net.nodes, net.dtype);
            epsilon    = net.epsilon;
            layer      = net.layer;
            lHL        = numLayers-1;
            numParamsPerlayer_2 = net.numParamsPerlayer_2;
            smlIdx     = net.similarIdx;
            
            deltaM     = cell(numLayers, 1);
            deltaS     = cell(numLayers, 1);
            deltaMx    = cell(numLayers, 1);
            deltaSx    = cell(numLayers, 1);
            deltaMxs   = cell(numLayers, 1);
            deltaMdxs  = cell(numLayers, 1);
            deltaSxs   = cell(numLayers, 1);
            deltaSdxs  = cell(numLayers, 1);                      
            if net.lastLayerUpdate
                if net.learnSv == 0
                    % Update hidden states for the last hidden layer
                    if isempty(Sy)
                        R = net.sv.^2;
                    else
                        R = net.sv.^2; %+ Sy;
                    end                   
                    if isempty(udIdx)
                        Szv = Sa{end} + R;
                        if net.imperfect_obs
                            deltaMz = y - ma{lHL+1};
                            deltaSz = Sy - Szv;
                            % Innovation vector
                            [deltaMz, deltaSz] = tagi.inovationVector(Szv,...
                                deltaMz, deltaSz, net.gpu);
                            % RTS update
                            deltaMz = Sa{end} .* deltaMz;
                            deltaSz = Sa{end} .* deltaSz .* Sa{end};
                            
%                             [temp_deltaMz,...
%                              temp_deltaSz] = tagi.fowardHiddenStateUpdate(ma{lHL+1},...
%                              Szv, J{lHL+1}.*Sz{lHL+1}, y, net.gpu);
                        else
                            [deltaMz,...
                             deltaSz] = tagi.fowardHiddenStateUpdate(ma{lHL+1},...
                             Szv, J{lHL+1}.*Sz{lHL+1}, y, net.gpu);
                        end
                    else
                        mzf = ma{end}(udIdx);
                        Szf = J{lHL+1}(udIdx) .* Sz{lHL+1}(udIdx);
                        ys  = y;
                        Szv = Sa{end}(udIdx) + R;
                        deltaMz = zeros(size(mz{lHL+1}), 'like', mz{lHL+1});
                        deltaSz = zeros(size(Sz{lHL+1}), 'like', Sz{lHL+1});
                        [deltaMz(udIdx),...
                         deltaSz(udIdx)] = tagi.fowardHiddenStateUpdate(mzf,...
                         Szv, Szf, ys, net.gpu);
                    end
                elseif net.learnSv==1                   
                    if strcmp(net.task, 'regression') ...
                            && strcmp(net.noiseType, 'hete')
                        [mla, mv2a] = tagi.detachMeanVar(ma{end}, net.nl,...
                            net.nv2, B, rB);
                        [Sla, Sv2a] = tagi.detachMeanVar(Sa{end}, net.nl,...
                            net.nv2, B, rB);
                        [Slz, ~]  = tagi.detachMeanVar(Sz{end}, net.nl,...
                            net.nv2, B, rB);
                        [Jl, Jv2] = tagi.detachMeanVar(J{end}, net.nl,...
                            net.nv2, B, rB);
                        % Activate log(\sigma_v2)
                        [mv2a, Sv2a, Cv2a] = act.expFun(mv2a, Sv2a, net.gpu);
                        
                        [deltaMlz, deltaSlz, deltaMv2z,...
                            deltaSv2z] = tagi.noiseUpdate4regression(Slz,...
                            mla, Sla, Jl, Jv2, mv2a, Sv2a, Cv2a, y, Sy, net.sv,...
                            net);
                        deltaMz = tagi.attachMeanVar(deltaMlz, deltaMv2z,...
                            net.nl, net.nv2, B, rB);
                        deltaSz = tagi.attachMeanVar(deltaSlz, deltaSv2z,...
                            net.nl, net.nv2, B, rB);
                    elseif strcmp(net.task, 'regression') ...
                            && strcmp(net.noiseType, 'full')
                        [mla, mLa] = tagi.detachMeanVar(ma{end}, net.nl,...
                            net.nLchol, B, rB);
                        [Sla, SLa] = tagi.detachMeanVar(Sa{end}, net.nl,...
                            net.nLchol, B, rB);
                        % transform the diagonal elements into positive domain
                        [mLa_, SLa_, CLa_] = agvi.transform_chol_vec(mLa, SLa, net.gpu);
                        % retrieve post parameters for l and z
                        [deltaMlz, deltaSlz, mLa_post_, SLa_post_] = agvi.full_noiseUpdate4regression(mla, Sla, mLa_, SLa_, y);
                        % transform L back to the real domain
                        [mLa_post, SLa_post] = agvi.full_noiseBackwardUpdate(mLa, SLa, mLa_, SLa_, CLa_, mLa_post_, SLa_post_, net.gpu);

                        deltaMLz = mLa_post - mLa;
                        deltaSLz = SLa_post - SLa;

                        deltaMz = tagi.attachMeanVar(deltaMlz, deltaMLz,...
                            net.nl, net.nLchol, B, rB);
                        deltaSz = tagi.attachMeanVar(deltaSlz, deltaSLz,...
                            net.nl, net.nLchol, B, rB);
                    elseif strcmp(net.task, 'regression') ...
                            && strcmp(net.noiseType, 'homo')
                        mv2a = net.sv(1);
                        Sv2a = net.sv(2);
                        mla  = ma{end};
                        Slz  = Sz{end};
                        Sla  = Sa{end};
                        Jl   = J{end};  
                        [deltaMz, deltaSz, deltaMv2z,...
                         deltaSv2z] = tagi.homoNoiseUpdate4regression(Slz,...
                         mla, Sla, Jl, mv2a, Sv2a, y, net.gpu);
                        net.sv(1) = net.sv(1) + sum(deltaMv2z, 1);
                        net.sv(2) = net.sv(2) + sum(deltaSv2z, 1);                       
                    end                    
                end
            else
                deltaMz = y;
                deltaSz = Sy;
            end
            sv = net.sv;
            for k = (numLayers-1) : -1 : 1
                if kernelSize(k) == stride(k) ...
                        || (kernelSize(k) == imgW(k) ...
                        && stride(k)==1)
                    overlap = 0;
                else
                    overlap = 1; 
                end
                if isempty(mdxs{k+1}); nSz = Sz{k+1}; else; nSz = Sdxs{k+1}; end
                if isempty(mdxs{k}); cSz = Sz{k}; else; cSz = Sdxs{k}; end
                
                cSxs = Sxs{k};
                idxw = (numParamsPerlayer_2(1, k)+1):numParamsPerlayer_2(1, k+1);
                %Shortcut connection for residual network
                if net.xsc(k+1)~=0 ...
                        && (net.filter(net.xsc(k+1)) ~= net.filter(k+1) ...
                        || net.imgW(net.xsc(k+1)) ~= net.imgH(k+1))
                    [deltaMx{k+1},...
                        deltaSx{k+1}] = tagi.inovationVector(Sxs{k+1},...
                        deltaMzx, deltaSzx, net.gpu);
                    idxXsc = net.xsc(k+1);  
                    idxwx = (numParamsPerlayer_2(3, idxXsc)+1):numParamsPerlayer_2(3, idxXsc+1);
                    if idxXsc>1                                 
                        [deltaMxs{idxXsc}, deltaSxs{idxXsc}, deltaMdxs{idxXsc},...
                         deltaSdxs{idxXsc}] = tagi.xshortDelta(deltaMx{k+1},...
                         deltaSx{k+1}, Sxs{idxXsc}, Sdxs{idxXsc}, J{idxXsc},...
                         mwx(idxwx), net.idxSzzUdXsc{idxXsc},...
                         net.idxFCzwaXsc(idxXsc, :), filter(idxXsc), B, rB,...
                         size(net.idxFCzwaXsc{idxXsc, 2}, 1), net.gpu); 
                    end                   
                elseif net.xsc(k+1)~=0 ...
                        && (net.filter(net.xsc(k+1)) == net.filter(k+1) ...
                        ||net.imgW(net.xsc(k+1)) == net.imgH(k+1))
                    [deltaMx{k+1},...
                        deltaSx{k+1}] = tagi.inovationVector(Sxs{k+1},...
                        deltaMzx, deltaSzx, net.gpu);
                    idxXsc = net.xsc(k+1);
                    if idxXsc > 1 && ~isempty(Sxs{idxXsc})                      
                        [deltaMxs{idxXsc}, deltaSxs{idxXsc}, deltaMdxs{idxXsc},...
                         deltaSdxs{idxXsc}] = tagi.xshortDelta(deltaMx{k+1},...
                         deltaSx{k+1}, Sxs{idxXsc}, Sdxs{idxXsc}, J{idxXsc},...
                         [],  [], [], [], [], rB, [], net.gpu);
                    elseif idxXsc > 1 ...
                            && isempty(Sdxs{idxXsc}) ...
                            && isempty(Sxs{idxXsc}) % First shortcut
                        [~, ~, deltaMdxs{idxXsc},...
                         deltaSdxs{idxXsc}] = tagi.xshortDelta(deltaMx{k+1},...
                         deltaSx{k+1}, [], Sz{idxXsc}, J{idxXsc}, [], [],...
                         [], [], [], rB, [], net.gpu);
                    end
                end   
                
                % Innovation vector
                [deltaM{k+1}, deltaS{k+1}] = tagi.inovationVector(nSz,...
                    deltaMz, deltaSz, net.gpu);
                
                % Max pooling 
                if layer(k+1) == net.layerEncoder.mp       
                    [deltaMz, deltaSz, deltaMzx,...
                        deltaSzx] = tagi.mpHiddenStateBackwardPass(cSz, cSxs,...
                        J{k}, deltaM{k+1}, deltaS{k+1}, maxIdx{k+1}, rB,...
                        overlap, net.gpu);
                    
                % Average pooling     
                elseif layer(k+1) == net.layerEncoder.ap 
                    [deltaMz, deltaSz, deltaMzx,...
                        deltaSzx] = tagi.agHiddenStateBackwardPass(cSz, cSxs,...
                        J{k}, size(net.idxPooling{k}, 2), deltaM{k+1},...
                        deltaS{k+1}, net.idxSzzUd{k}, imgW(k+1), imgH(k+1),...
                        filter(k+1), kernelSize(k), B, rB, overlap, net.gpu);
                    
                % Convolutional     
                elseif layer(k+1) == net.layerEncoder.conv 
                    if k > 1||net.convariateEstm
                        if B == 1 && rB == 1
                            [deltaMz, deltaSz, deltaMzx,...
                             deltaSzx] = tagi.convHiddenStateBackwardPassB1(cSz,...
                             cSxs, J{k}, mw(idxw), deltaM{k+1}, deltaS{k+1},...
                            net.idxSzzUd{smlIdx(k)},...
                            net.idxFCzwa(smlIdx(k), :),...
                            imgW(k), imgH(k), filter(k), net.gpu);
                        else
                            [deltaMz, deltaSz, deltaMzx,...
                             deltaSzx] = tagi.convHiddenStateBackwardPass(cSz,...
                             cSxs, J{k}, mw(idxw), deltaM{k+1}, deltaS{k+1},...
                             net.idxSzzUd{smlIdx(k)},...
                             net.idxFCzwa(smlIdx(k), :), imgW(k), imgH(k),...
                             filter(k), B, rB, net.gpu);
                        end                       
                    end
                    
                % Transposed convolutional
                elseif layer(k+1) == net.layerEncoder.tconv 
                    if k > 1 || net.convariateEstm
                        [deltaMz, deltaSz, deltaMzx,...
                         deltaSzx] = tagi.tconvHiddenStateBackwardPass(cSz,...
                         cSxs, J{k}, mw(idxw), deltaM{k+1}, deltaS{k+1},...
                         net.idxSzzUd{k}, net.idxFCzwa(k, :), imgW(k),...
                         imgH(k), filter(k), B, rB, net.gpu);                       
                    end
                    
                % Normalization     
                elseif layer(k+1) == net.layerEncoder.ln ...
                        || layer(k+1) == net.layerEncoder.bn                     
                    if k > 1 || net.convariateEstm
                        Shat = tagi.distributeNormMeanVar(Sra{k}, nodes(k),...
                            imgW(k), imgH(k), filter(k), B, rB, layer(k),...
                            layer(k+1), net.layerEncoder);
                        [deltaMz, deltaSz, deltaMzx,...
                            deltaSzx] = tagi.normHiddenStateBackwardPass(cSz,...
                            cSxs, J{k}, mw(idxw), Shat, epsilon, deltaM{k+1},...
                            deltaS{k+1}, imgW(k), imgH(k), filter(k), B, rB,...
                            layer(k), net.layerEncoder, net.gpu);
                    end 
                    
                % Full-connected     
                elseif  layer(k+1) == net.layerEncoder.fc
                    if k > 1||net.convariateEstm
                        [deltaMz, deltaSz, deltaMzx,...
                         deltaSzx] = tagi.fcHiddenStateBackwardPass(cSz,...
                         cSxs, J{k}, mw(idxw), deltaM{k+1}, deltaS{k+1},...
                         nodes(k), nodes(k+1), B, rB, net.gpu);
                    end                  
                end
                
                % Update hidden states from shortcut
                if ~isempty(deltaMxs{k}) && ~isempty(deltaMdxs{k})
                    [deltaMzx, deltaSzx, deltaMz,...
                        deltaSz] = arrayfun(@fourPlus, deltaMzx, deltaSzx,...
                        deltaMz, deltaSz, deltaMxs{k}, deltaSxs{k},...
                        deltaMdxs{k}, deltaSdxs{k});
                elseif ~isempty(deltaMdxs{k}) && isempty(deltaMxs{k})
                    [deltaMz, deltaSz] = arrayfun(@twoPlus, deltaMz, deltaSz,...
                        deltaMdxs{k}, deltaSdxs{k});
                end
            end
            deltaMz0 = deltaMz;
            deltaSz0 = deltaSz;        
        end

        function deltaTheta = parameterBackwardPass(net, theta, normStat,...
                states, deltaM, deltaS, deltaMx, deltaSx)
            % Initialization
            [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = tagi.extractParameters(theta);
            [~, ~, ma, ~, ~, ~, ~, ~, ~] = tagi.extractStates(states);
            [mra, Sra] = tagi.extractNormStat(normStat);
            numLayers  = length(net.nodes);
            imgW       = cast(net.imgW, net.dtype);
            imgH       = cast(net.imgH, net.dtype);
            filter     = cast(net.filter, net.dtype);
            kernelSize = cast(net.kernelSize, net.dtype);
            B          = cast(net.batchSize, net.dtype);
            rB         = cast(net.repBatchSize, net.dtype);
            nodes      = cast(net.nodes, net.dtype);
            epsilon    = net.epsilon;
            layer      = net.layer;
            numParamsPerlayer_2 = net.numParamsPerlayer_2;
            smlIdx     = net.similarIdx;
            
            deltaMw    = mw;
            deltaSw    = Sw;
            deltaMb    = mb;
            deltaSb    = Sb;
            deltaMwx   = mwx;
            deltaSwx   = Swx;
            deltaMbx   = mbx;
            deltaSbx   = Sbx;
            for k = (numLayers - 1) : -1 : 1
                idxw = (numParamsPerlayer_2(1, k)+1):numParamsPerlayer_2(1, k+1);
                idxb = (numParamsPerlayer_2(2, k)+1):numParamsPerlayer_2(2, k+1);
                %Shortcut connection for residual network
                if net.xsc(k+1)~=0 ...
                        && (net.filter(net.xsc(k+1)) ~= net.filter(k+1) ...
                        ||net.imgW(net.xsc(k+1)) ~= net.imgH(k+1))
                    idxXsc = net.xsc(k+1); 
                    idxwx = (numParamsPerlayer_2(3, idxXsc)+1):numParamsPerlayer_2(3, idxXsc+1);
                    idxbx = (numParamsPerlayer_2(4, idxXsc)+1):numParamsPerlayer_2(4, idxXsc+1);
                    
                    [deltaMwx(idxwx), deltaSwx(idxwx), deltaMbx(idxbx),...
                     deltaSbx(idxbx)] = tagi.convParameterBackwardPass(deltaMwx(idxwx),...
                     deltaSwx(idxwx), deltaMbx(idxbx), deltaSbx(idxbx),...
                     Swx(idxwx), Sbx(idxbx), ma{idxXsc}, deltaMx{k+1},...
                     deltaSx{k+1}, net.idxFmwaXsc(idxXsc, :),...
                     net.paddingXsc(idxXsc), 1, filter(idxXsc), imgW(k+1),...
                     imgH(k+1), filter(k+1), B, rB, net.gpu);
                end 
                
                % Convolutional     
                if layer(k+1) == net.layerEncoder.conv  
                    if B == 1 && rB == 1
                        [deltaMw(idxw), deltaSw(idxw), deltaMb(idxb),...
                         deltaSb(idxb)] = tagi.convParameterBackwardPassB1(Sw(idxw),...
                         Sb(idxb), ma{k}, deltaM{k+1}, deltaS{k+1},...
                         net.idxFmwa(smlIdx(k), :), net.padding(k),...
                         kernelSize(k), filter(k), imgW(k+1), imgH(k+1),...
                         filter(k+1), net.gpu);
                    else
                        [deltaMw(idxw), deltaSw(idxw), deltaMb(idxb),...
                         deltaSb(idxb)] = tagi.convParameterBackwardPass(deltaMw(idxw),...
                         deltaSw(idxw), deltaMb(idxb), deltaSb(idxb),...
                         Sw(idxw), Sb(idxb), ma{k}, deltaM{k+1}, deltaS{k+1},...
                         net.idxFmwa(smlIdx(k), :), net.padding(k),...
                         kernelSize(k), filter(k), imgW(k+1), imgH(k+1),...
                         filter(k+1), B, rB, net.gpu);
                    end                                                         
                % Transposed convolutional
                elseif layer(k+1) == net.layerEncoder.tconv 
                    [deltaMw(idxw), deltaSw(idxw), deltaMb(idxb),...
                     deltaSb(idxb)] = tagi.tconvParameterBackwardPass(deltaMw(idxw),...
                     deltaSw(idxw), deltaMb(idxb), deltaSb(idxb), Sw(idxw),...
                     Sb(idxb), ma{k}, deltaM{k+1}, deltaS{k+1},...
                     net.idxSwzUd{k}, net.idxFCwz(k, :), kernelSize(k),...
                     filter(k), imgW(k+1), imgH(k+1), filter(k+1), B, rB,...
                     net.gpu); 
                    
                % Normalization     
                elseif layer(k+1) == net.layerEncoder.ln ...
                        || layer(k+1) == net.layerEncoder.bn  
                    mhat = tagi.distributeNormMeanVar(mra{k}, nodes(k),...
                        imgW(k), imgH(k), filter(k), B, rB, layer(k), layer(k+1), net.layerEncoder);
                    Shat = tagi.distributeNormMeanVar(Sra{k}, nodes(k),...
                        imgW(k), imgH(k), filter(k), B, rB, layer(k),...
                        layer(k+1), net.layerEncoder);
                    [deltaMw(idxw), deltaSw(idxw), deltaMb(idxb),...
                     deltaSb(idxb)] = tagi.normParameterBackwardPass(deltaMw(idxw),...
                     deltaSw(idxw), deltaMb(idxb), deltaSb(idxb), Sw(idxw),...
                     Sb(idxb), ma{k}, mhat, Shat, epsilon, deltaM{k+1},...
                     deltaS{k+1}, nodes(k), imgW(k), imgH(k), filter(k), B,...
                     rB, layer(k), net.layerEncoder, net.gpu);        
                    
                % Full-connected     
                elseif  layer(k+1) == net.layerEncoder.fc 
                    [deltaMw(idxw), deltaSw(idxw), deltaMb(idxb),...
                     deltaSb(idxb)] = tagi.fcParameterBackwardPass(deltaMw(idxw),...
                     deltaSw(idxw), deltaMb(idxb), deltaSb(idxb), Sw(idxw),...
                     Sb(idxb), ma{k}, deltaM{k+1}, deltaS{k+1}, nodes(k),...
                     nodes(k+1), B, rB, net.gpu);
                end
            end
            deltaTheta = tagi.compressParameters(deltaMw, deltaSw, deltaMb,...
                deltaSb, deltaMwx, deltaSwx, deltaMbx, deltaSbx);           
        end             
        
        % Pooling layer
        function [mz, Sz, maxIdx] = mpMeanVar(mz, Sz, maS, ma, Sa,...
                idxpooling, maxIdx, rB, gpu)
            n = size(idxpooling, 1);   
            for t = 1:rB
                maSloop = maS(:, t);
                [~, idx] = max(maSloop(idxpooling), [], 2);
                if gpu
                    n   = gpuArray(cast(n, 'int32'));
                    col = gpuArray.colon(1, n);
                    col = col(:);
                    fun = @(x, y, z) (x-1).*y + z;
                    idx = arrayfun(fun, idx, n, col);
                else
                    col = colon(1,n)';
                    idx = (idx-1)*n + col;
                end
                maxIdx(:, t) = idxpooling(idx);
                mz(:, t) = ma(maxIdx(:, t), t);
                Sz(:, t) = Sa(maxIdx(:, t), t);
            end
        end
        function [mz, Sz] = apMeanVar(mz, Sz, ma, Sa, idxPooling, padding,...
                rB)
            n   = size(idxPooling, 2);
            if padding ~= 0
                zeroPad = zeros(1, size(ma, 2), 'like', ma);
                ma = [ma; zeroPad];
                Sa = [Sa; zeroPad];
            end
            for t = 1:rB  
                maloop = ma(:, t);
                Saloop = Sa(:, t);
                mz(:, t) = mean(maloop(idxPooling), 2);
                Sz(:, t) = sum(Saloop(idxPooling), 2)./(n^2);
            end           
        end             
        function [deltaMz, deltaSz,...
                deltaMxs, deltaSxs] = mpHiddenStateBackwardPass(Sz, Sxs, J,...
                deltaM, deltaS, maxIdx, rB, overlap, gpu)
            deltaMz  = Sz;
            deltaSz  = Sz;
            deltaMxs = Sxs;
            deltaSxs = Sxs;
            n = single(size(Sz, 1));
            if gpu
                if isempty(Sxs)
                    for t = 1:rB
                        Czz = bsxfun(@times, J(:, t), Sz(:, t));
                        Czz = Czz(maxIdx(:, t));
                        if overlap == 1
                            [deltaMzloop,...
                             deltaSzloop] = arrayfun(@vectorizedDelta, Czz,...
                             deltaM(:, t), deltaS(:, t));
                            deltaMz(:, t) = accumarray(maxIdx(:, t),...
                             deltaMzloop, [n, 1], @sum);
                            deltaSz(:, t) = accumarray(maxIdx(:, t),...
                             deltaSzloop , [n, 1], @sum);
                        else
                            [deltaMz(maxIdx(:, t), t),...
                              deltaSz(maxIdx(:, t), t)] = arrayfun(@vectorizedDelta,...
                              Czz, deltaM(:, t), deltaS(:, t));
                        end
                    end
                else
                    for t = 1:rB
                        if overlap == 1
                            [deltaMzloop, deltaSzloop, deltaMxsloop,...
                              deltaSxsloop] = arrayfun(@vectorized4delta,...
                              J(maxIdx(:, t), t), Sz(maxIdx(:, t), t),...
                              Sxs(maxIdx(:, t), t), deltaM(:, t), deltaS(:, t));
                            deltaMz(:, t)  = accumarray(maxIdx(:, t),...
                                deltaMzloop, [n, 1], @sum);
                            deltaSz(:, t)  = accumarray(maxIdx(:, t),...
                                deltaSzloop , [n, 1], @sum);
                            deltaMxs(:, t) = accumarray(maxIdx(:, t),...
                                deltaMxsloop, [n, 1], @sum);
                            deltaSxs(:, t) = accumarray(maxIdx(:, t),...
                                deltaSxsloop , [n, 1], @sum);
                        else
                            [deltaMz(maxIdx(:, t), t),...
                              deltaSz(maxIdx(:, t), t),...
                              deltaMxs(maxIdx(:, t), t),...
                              deltaSxs(maxIdx(:, t), t)] = arrayfun(@vectorized4delta,...
                              J(maxIdx(:, t), t), Sz(maxIdx(:, t), t),...
                              Sxs(maxIdx(:, t), t), deltaM(:, t), deltaS(:, t));
                        end
                    end
                end
            else
                if isempty(Sxs)
                    for t = 1:rB
                        Czz = J(:, t).*Sz(:, t);
                        Czz = Czz(maxIdx(:, t));
                        if overlap == 1
                            deltaMzloop   = Czz.*deltaM(:, t);
                            deltaSzloop   = Czz.*deltaS(:, t) .* Czz;
                            deltaMz(:, t) = accumarray(maxIdx(:, t), ...
                                deltaMzloop, [n, 1], @sum);
                            deltaSz(:, t) = accumarray(maxIdx(:, t), ...
                                deltaSzloop , [n, 1], @sum);
                        else
                            deltaMz(maxIdx(:, t), t) = Czz .* deltaM(:, t);
                            deltaSz(maxIdx(:, t), t) = Czz .* deltaS(:, t) .* Czz;
                        end
                    end
                else
                    for t = 1:rB
                        Czz = J(:, t) .* Sz(:, t);
                        Czz = Czz(maxIdx(:, t));
                        Czx = J(:, t) .* Sxs(:, t);
                        Czx = Czx(maxIdx(:, t));
                        if overlap == 1
                            deltaMzloop    = Czz .* deltaM(:, t);
                            deltaSzloop    = Czz .* deltaS(:, t) .* Czz;
                            deltaMxsloop   = Czx .* deltaM(:, t);
                            deltaSxsloop   = Czx .* deltaS(:, t) .* Czx;
                            deltaMz(:, t)  = accumarray(maxIdx(:, t),...
                                deltaMzloop, [n, 1], @sum);
                            deltaSz(:, t)  = accumarray(maxIdx(:, t),...
                                deltaSzloop, [n, 1], @sum);
                            deltaMxs(:, t) = accumarray(maxIdx(:, t),...
                                deltaMxsloop, [n, 1], @sum);
                            deltaSxs(:, t) = accumarray(maxIdx(:, t),...
                                deltaSxsloop, [n, 1], @sum);
                        else
                            deltaMz(maxIdx(:, t), t) = Czz .* deltaM(:, t);
                            deltaSz(maxIdx(:, t), t) = Czz .* deltaS(:, t) .* Czz;
                        end
                    end
                end
            end
        end              
        function [deltaMz, deltaSz,...
                deltaMxs, deltaSxs] = agHiddenStateBackwardPass(Sz, Sxs, J,...
                n, deltaM, deltaS, idx,  wo, ho, fo, ki, B, rB, overlap, gpu)    
            deltaMz  = Sz;
            deltaSz  = Sz;
            deltaMxs = Sxs;
            deltaSxs = Sxs;
            n = cast(n, 'like', Sz);
            if gpu
                if isempty(Sxs)
                    for t = 1:rB
                        if overlap == 0
                            deltaMzloop = reshape(repmat(reshape(repmat(transpose(deltaM(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                            deltaSzloop = reshape(repmat(reshape(repmat(transpose(deltaS(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                        else
                            zeroPadding = zeros(1,1,'like',deltaM);
                            deltaMzloop = [deltaM(:, t); zeroPadding];
                            deltaSzloop = [deltaS(:, t); zeroPadding];
                            deltaMzloop = deltaMzloop(idx);
                            deltaSzloop = deltaSzloop(idx);
                        end
                        [deltaMzloop,...
                          deltaSzloop] = arrayfun(@vectorizedDelta_V2,...
                          J(:, t), Sz(:, t)/n, deltaMzloop, deltaSzloop);
                        deltaMz(:, t) = sum(deltaMzloop, 2);
                        deltaSz(:, t) = sum(deltaSzloop, 2);
                    end
                else
                    Czz = bsxfun(@times, J, Sz);
                    Czx = bsxfun(@times, J, Sxs);
                    for t = 1:rB
                        if overlap == 0
                            deltaMloop = reshape(repmat(reshape(repmat(transpose(deltaM(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                            deltaSloop = reshape(repmat(reshape(repmat(transpose(deltaS(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo *B ]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                        else
                            zeroPadding = zeros(1,1,'like',deltaM);
                            deltaMloop = [deltaM(:, t); zeroPadding];
                            deltaSloop = [deltaS(:, t); zeroPadding];
                            deltaMloop = deltaMloop(idx);
                            deltaSloop = deltaSloop(idx);
                        end
                        [deltaMzloop, deltaSzloop, deltaMxsloop,...
                         deltaSxsloop] = arrayfun(@vectorized4delta, 1/n,...
                         Czz(:, t), Czx(:, t), deltaMloop, deltaSloop);
                        deltaMz(:, t) = sum(deltaMzloop, 2);
                        deltaSz(:, t) = sum(deltaSzloop, 2);
                        deltaMxs(:, t) = sum(deltaMxsloop, 2);
                        deltaSxs(:, t) = sum(deltaSxsloop, 2);
                    end
                end
            else
                if isempty(Sxs)
                    for t = 1:rB
                        Czz = (J(:, t).*Sz(:, t))/n;
                        if overlap == 0
                            deltaMzloop = reshape(repmat(reshape(repmat(transpose(deltaM(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                            deltaSzloop = reshape(repmat(reshape(repmat(transpose(deltaS(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                        else
                            zeroPadding = zeros(1,1,'like',deltaM);
                            deltaMzloop = [deltaM(:, t); zeroPadding];
                            deltaSzloop = [deltaS(:, t); zeroPadding];
                            deltaMzloop = deltaMzloop(idx);
                            deltaSzloop = deltaSzloop(idx);
                        end
                        deltaMzloop   = Czz .* deltaMzloop;
                        deltaSzloop   = Czz .* deltaSzloop .* Czz;
                        deltaMz(:, t) = sum(deltaMzloop, 2);
                        deltaSz(:, t) = sum(deltaSzloop, 2);
                    end
                else
                    Czz = (J.*Sz)/n;
                    Czx = (J.*Sxs)/n;
                    for t = 1:rB                       
                        if overlap == 0
                            deltaMloop = reshape(repmat(reshape(repmat(transpose(deltaM(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                            deltaSloop = reshape(repmat(reshape(repmat(transpose(deltaS(:, t)),...
                                [ki, 1]), [ki * ho, wo * fo * B]), [ki, 1]),...
                                [ho * wo * fo * ki * ki * B, 1]);
                        else
                            zeroPadding = zeros(1,1,'like',deltaM);
                            deltaMloop = [deltaM(:, t); zeroPadding];
                            deltaSloop = [deltaS(:, t); zeroPadding];
                            deltaMloop = deltaMloop(idx);
                            deltaSloop = deltaSloop(idx);
                        end
                        deltaMzloop    = Czz(:, t) .* deltaMloop;
                        deltaSzloop    = Czz(:, t) .* deltaSloop .* Czz(:, t);
                        deltaMxsloop   = Czx(:, t) .* deltaMloop;
                        deltaSxsloop   = Czx(:, t) .* deltaSloop .* Czx(:, t);
                        deltaMz(:, t)  = sum(deltaMzloop, 2);
                        deltaSz(:, t)  = sum(deltaSzloop, 2);
                        deltaMxs(:, t) = sum(deltaMxsloop, 2);
                        deltaSxs(:, t) = sum(deltaSxsloop, 2);
                    end
                end
            end
        end        
        function [mz, Sz, maxIdx] = mpMeanVarMexcuda(maS, ma, Sa,...
                idxPooling, wo, ho, fo, wi, hi, fi, ki, B)
            [mz, Sz, maxIdx] = mpHiddenStateForwardPass4matlab(maS, ma, Sa,...
                idxPooling, wo, ho, fo, wi, hi, fi, ki, B, overlap); 
        end
        function [mz, Sz] = apMeanVarMexcuda(ma, Sa, idxPooling, wo, ho, fo,...
                wi, hi, fi, ki, B, overlap)
            [mz, Sz] = apHiddenStateForwardPass4matlab(ma, Sa, idxPooling,...
                wo, ho, fo, wi, hi, fi, ki, B, overlap);            
        end 
        function [deltaMz, deltaSz,...
                deltaMxs, deltaSxs] = apHiddenStateBackwardPassMexcuda(Sz,...
                Sxs, J, deltaM, deltaS, idxSzzUd, wo, ho, fo, wi, hi, fi,...
                ki, B, overlap) 
            [deltaMz, deltaSz] = apHiddenStateBackwardPass4matlab(Sz, J,...
             deltaM, deltaS, idxSzzUd, wo, ho, fo, wi, hi, fi, ki, B, overlap);  
            if ~isempty(Sxs)
                [deltaMxs, deltaSxs] = apHiddenStateBackwardPass4matlab(Sxs,...
                 J, deltaM, deltaS, idxSzzUd, wo, ho, fo, wi, hi, fi, ki, B,...
                 overlap);                
            else
                deltaMxs = [];
                deltaSxs = [];
            end
        end
        
        % Full connected layer 
        function [mz, Sz] = fcMeanVar(mz, Sz, mw, Sw, mb, Sb, ma, Sa, ni,...
                no, B, rB, gpu)
            idxSum = 1;
            if any(isnan(mb))
                mb = zeros(1,1,'like', mw);
                Sb = zeros(1,1,'like', Sw);               
            else
                mb = repmat(mb, [B, 1]);
                Sb = repmat(Sb, [B, 1]);
            end
            mw  = repmat(reshape(mw, [ni, no]), [1, B]);                     
            Sw  = repmat(reshape(Sw, [ni, no]), [1, B]);                                  
            if gpu
                for t = 1:rB
                    maloop = reshape(repmat(reshape(ma(:, t), [ni, B]),...
                        [no, 1]), [ni, no*B]);
                    Saloop = reshape(repmat(reshape(Sa(:, t), [ni, B]),...
                        [no, 1]), [ni, no*B]);
                    [mzloop, Szloop] = arrayfun(@vectorizedMeanVar, maloop,...
                        mw, Saloop, Sw);
                    mzloop = transpose(sum(mzloop, idxSum));
                    Szloop = transpose(sum(Szloop, idxSum));
                    [mz(:, t), Sz(:, t)] = arrayfun(@twoPlus, mzloop,...
                        Szloop, mb, Sb);
                end
            else
                for t = 1:rB
                    maloop = reshape(repmat(reshape(ma(:, t), [ni, B]),...
                        [no, 1]), [ni, no*B]);
                    Saloop = reshape(repmat(reshape(Sa(:, t), [ni, B]),...
                        [no, 1]), [ni, no*B]);
                    [mzloop, Szloop] = vectorizedMeanVar(maloop, mw, Saloop,...
                        Sw);
                    mzloop = transpose(sum(mzloop, idxSum));
                    Szloop = transpose(sum(Szloop, idxSum));
                    mz(:, t) = mzloop + mb;
                    Sz(:, t) = Szloop + Sb;
                end
            end            
        end
        function [mz, Sz] = fcMeanVarB1(mw, Sw, mb, Sb, ma, Sa, ni, no,...
                gpu)
            if any(isnan(mb))
                mb = zeros(1,1,'like', mw);
                Sb = zeros(1,1,'like', Sw);               
            end
            mw = reshape(mw, [ni, no]);                     
            Sw = reshape(Sw, [ni, no]); 
            if gpu
                [mzloop, Szloop] = arrayfun(@vectorizedMeanVar, ma, mw,...
                    Sa, Sw);
                mzloop = sum(mzloop, 1);
                Szloop = sum(Szloop, 1);
                mzloop = mzloop(:);
                Szloop = Szloop(:);
                [mz, Sz] = arrayfun(@twoPlus, mzloop, Szloop, mb, Sb);
            else
                [mzloop, Szloop] = vectorizedMeanVar(ma, mw, Sa, Sw);
                mzloop = transpose(sum(mzloop, 1));
                Szloop = transpose(sum(Szloop, 1));
                mz = mzloop + mb;
                Sz = Szloop + Sb;
            end            
        end  
        function [mz, Sz, Szf] = fcMeanCov(mw, Sw, mb, Sb, ma, Sa, Saf, ni,...
                no, B, gpu)
            idxSum = 1;
            if any(isnan(mb))
                mb = zeros(1,1,'like', mw);
                Sb = zeros(1,1,'like', Sw);               
            else
                mb = repmat(mb, [B, 1]);
                Sb = repmat(Sb, [B, 1]);
            end
            mwd = repmat(reshape(mw, [ni, no]), [1, B]);                     
            Swd = repmat(reshape(Sw, [ni, no]), [1, B]);  
            if gpu
                mad = reshape(repmat(reshape(ma, [ni, B]), [no, 1]),...
                    [ni, no * B]);
                Sad = reshape(repmat(reshape(Sa, [ni, B]), [no, 1]),...
                    [ni, no * B]);
                Sz  = Swd .* mad .* mad + Sad .* Swd;
                mz  = mad .* mwd;
                mz  = transpose(sum(mz, idxSum)) + mb;
                Sz  = transpose(sum(Sz, idxSum)) + Sb;
                
                mw1  = repmat(reshape(mw, [ni, 1, no]), [1, 1, no * B]);
                mw2  = repmat(reshape(repmat(reshape(mw, [ni, no]),...
                    [no, 1]), [1, ni, no * no]), [1 1 B]);
                mw3  = pagefun(@mtimes, mw1, mw2);
                Saf  = reshape(repmat(reshape(Saf, [ni, ni, B]),...
                    [1, no * no, 1]), [ni, ni, no * no * B]);
                Szf  = sum(sum(bsxfun(@times, mw3, Saf), 2), 1);
                Szf  = reshape(Szf, no * no, B);
                Szfd = Szf(1 : no + 1 : no * no, :);
                Sz   = Sz + Szfd(:);
                Szf(1 : no + 1 : no * no, :) = reshape(Sz, [no, B]); 
                Szf = Szf(:);
            else
                mad = reshape(repmat(reshape(ma, [ni, B]), [no, 1]),...
                    [ni, no * B]);
                Sad = reshape(repmat(reshape(Sa, [ni, B]), [no, 1]),...
                    [ni, no * B]);
                Sz  = Swd .* mad .* mad + Sad .* Swd;
                mz  = mad .* mwd;
                mz = transpose(sum(mz, idxSum)) + mb;
                Sz = transpose(sum(Sz, idxSum)) + Sb;
                
                mw1  = repmat(mw, [1, 1, no]);
                mw2  = repmat(reshape(mw, [1, ni, no]), [ni * no, 1, 1]);
                mw3  = mw2 .* mw1;
                mw3  = repmat(mw3, [1, 1, B]);
                Saf  = reshape(repmat(reshape(Saf, [ni, ni, B]),...
                    [no, no, 1]), [ni * no, ni, no * B]);
                Szf  = sum(reshape(sum(mw3 .* Saf, 2), [ni, 1, no * no * B]), 1);
                Szf  = reshape(Szf, no * no, B);
                Szfd = Szf(1 : no + 1 : no * no, :);
                Sz   = Sz + Szfd(:);
                Szf(1 : no + 1 : no * no, :) = reshape(Sz, [no, B]); 
                Szf = Szf(:);
            end
        end
        
        function [deltaMw, deltaSw,...
                deltaMb, deltaSb] = fcParameterBackwardPass(deltaMw,...
                deltaSw, deltaMb, deltaSb, Sw, Sb, ma, deltaMr, deltaSr,...
                ni, no, B, rB, gpu)  
            Cbz = repmat(Sb, [1, B]);
            if gpu 
                for t = 1:rB                   
                    maloop   = repmat(reshape(ma(:, t), [ni, B]), [no, 1]);               
                    deltaMrw = reshape(repmat(transpose(deltaMr(:, t)),...
                        [ni, 1]),[ni * no, B]);
                    deltaSrw = reshape(repmat(transpose(deltaSr(:, t)),...
                        [ni, 1]),[ni * no, B]);                  
                    % Weights
                    [deltaMrw, deltaSrw] = arrayfun(@vectorizedDelta_V2, Sw,...
                        maloop, deltaMrw, deltaSrw);
                    deltaMw(:, t) = sum(deltaMrw, 2);
                    deltaSw(:, t) = sum(deltaSrw, 2);
                    % Bias
                    if any(~isnan(Sb))                        
                        deltaMrb = reshape(deltaMr(:, t), [no, B]);
                        deltaSrb = reshape(deltaSr(:, t), [no, B]);                      
                        [deltaMrb, deltaSrb] = arrayfun(@vectorizedDelta,...
                            Cbz, deltaMrb, deltaSrb);
                        deltaMb(:, t) = sum(deltaMrb, 2);
                        deltaSb(:, t) = sum(deltaSrb, 2);
                    end
                end
            else
                for t = 1:rB
                    maloop   = repmat(reshape(ma(:, t), [ni, B]), [no, 1]);               
                    deltaMrw = reshape(repmat(transpose(deltaMr(:, t)),...
                        [ni, 1]),[ni*no, B]);
                    deltaSrw = reshape(repmat(transpose(deltaSr(:, t)),...
                        [ni, 1]),[ni*no, B]); 
                    Cwz      = Sw .* maloop;
                    deltaMrw = Cwz .* deltaMrw;
                    deltaSrw = Cwz .* deltaSrw .* Cwz;
                    deltaMw(:, t) = sum(deltaMrw, 2);
                    deltaSw(:, t) = sum(deltaSrw, 2);
                    if any(~isnan(Sb))
                        deltaMrb = reshape(deltaMr(:, t), [no, B]);
                        deltaSrb = reshape(deltaSr(:, t), [no, B]);                        
                        deltaMrb = Cbz .* deltaMrb;
                        deltaSrb = Cbz .* deltaSrb .* Cbz;
                        deltaMb(:, t) = sum(deltaMrb, 2);
                        deltaSb(:, t) = sum(deltaSrb, 2);
                    end
                end
            end  
            deltaMw = sum(deltaMw, 2);
            deltaSw = sum(deltaSw, 2);
            deltaMb = sum(deltaMb, 2);
            deltaSb = sum(deltaSb, 2);
        end
        function [deltaMw, deltaSw,...
                deltaMb, deltaSb] = fcParameterBackwardPassB1(Sw, Sb, ma,...
                deltaMr, deltaSr, ni, no, gpu)  
            Cbz      = Sb;                   
            maloop   = repmat(ma, [no, 1]);
            deltaMrw = repmat(transpose(deltaMr), [ni, 1]);
            deltaMrw = deltaMrw(:);
            deltaSrw = repmat(transpose(deltaSr), [ni, 1]);
            deltaSrw = deltaSrw(:);
            % Weights
            if gpu
                [deltaMrw, deltaSrw] = arrayfun(@vectorizedDelta_V2, Sw,...
                    maloop, deltaMrw, deltaSrw);
            else
                Cwa = Sw .* maloop;
                deltaMrw = Cwa .* deltaMrw;
                deltaSrw = (Cwa.^2) .* deltaSrw;                
            end
            deltaMw = sum(deltaMrw, 2);
            deltaSw = sum(deltaSrw, 2);
            % Bias
            if any(~isnan(Sb))
                if gpu
                    [deltaMrb, deltaSrb] = arrayfun(@vectorizedDelta, Cbz,...
                        deltaMr, deltaSr);
                else
                    [deltaMrb, deltaSrb] = vectorizedDelta(Cbz, deltaMr, deltaSr);
                end
                deltaMb = sum(deltaMrb, 2);
                deltaSb = sum(deltaSrb, 2);
            else
                deltaMb = Sb;
                deltaSb = Sb;
            end
        end
        function [deltaMz, deltaSz,...
                deltaMzx, deltaSzx] = fcHiddenStateBackwardPass(Sz, Sxs,...
                J, mw, deltaM, deltaS, ni, no, B, rB, gpu) 
            deltaMz  = Sz;
            deltaSz  = Sz;
            deltaMzx = Sxs;
            deltaSzx = Sxs;
            mw = repmat(reshape(mw, [ni, no]), [B, 1]);              
            if gpu
                Caz = bsxfun(@times, J, Sz);
                if isempty(Sxs)
                    for t = 1 : rB
                        deltaMzloop = reshape(repmat(reshape(deltaM(:, t),...
                            [no, B]), [ni, 1]), [no, ni * B])';
                        deltaSzloop = reshape(repmat(reshape(deltaS(:, t),...
                            [no, B]), [ni, 1]), [no, ni * B])';
                        [deltaMzloop,...
                         deltaSzloop] = arrayfun(@vectorizedDelta_V2, mw,...
                         Caz(:, t), deltaMzloop, deltaSzloop);
                        deltaMz(:, t) = sum(deltaMzloop, 2);
                        deltaSz(:, t) = sum(deltaSzloop, 2);
                    end
                else
                    Caxs = bsxfun(@times, J, Sxs);
                    for t = 1 : rB
                        deltaMloop = reshape(repmat(reshape(deltaM(:, t),...
                            [no, B]), [ni, 1]), [no, ni * B])';
                        deltaSloop = reshape(repmat(reshape(deltaS(:, t),...
                            [no, B]), [ni, 1]), [no, ni * B])';
                        [deltaMzloop, deltaSzloop, deltaMxsloop,...
                         deltaSxsloop] = arrayfun(@vectorized4delta, mw,...
                         Caz(:, t), Caxs(:, t), deltaMloop, deltaSloop);
                        deltaMz(:, t)  = sum(deltaMzloop, 2);
                        deltaSz(:, t)  = sum(deltaSzloop, 2);
                        deltaMzx(:, t) = sum(deltaMxsloop, 2);
                        deltaSzx(:, t) = sum(deltaSxsloop, 2);
                    end
                end
            else
                if isempty(Sxs)
                    for t = 1 : rB
                        Czz = J(:, t) .* Sz(:, t) .* mw;
                        deltaMzloop = reshape(repmat(reshape(deltaM(:, t),...
                            [no, B]), [ni, 1]), [no, ni*B])';
                        deltaSzloop = reshape(repmat(reshape(deltaS(:, t),...
                            [no, B]), [ni, 1]), [no, ni*B])';
                        deltaMzloop = Czz .* deltaMzloop;
                        deltaSzloop = Czz .* deltaSzloop .* Czz;
                        deltaMz(:, t) = sum(deltaMzloop, 2);
                        deltaSz(:, t) = sum(deltaSzloop, 2);
                    end
                else
                    for t = 1:rB
                        Czz = J(:, t) .*Sz (:, t) .* mw;
                        Czx = J(:, t) .*Sz (:, t) .* mw;
                        deltaMloop     = reshape(repmat(reshape(deltaM(:, t),...
                            [no, B]), [ni, 1]), [no, ni*B])';
                        deltaSloop     = reshape(repmat(reshape(deltaS(:, t),...
                            [no, B]), [ni, 1]), [no, ni*B])';
                        deltaMzloop    = Czz .* deltaMloop;
                        deltaSzloop    = Czz .* deltaSloop .* Czz;
                        deltaMxsloop   = Czx .* deltaMloop;
                        deltaSxsloop   = Czx .* deltaSloop .* Czx;
                        deltaMz(:, t)  = sum(deltaMzloop, 2);
                        deltaSz(:, t)  = sum(deltaSzloop, 2);
                        deltaMzx(:, t) = sum(deltaMxsloop, 2);
                        deltaSzx(:, t) = sum(deltaSxsloop, 2);
                    end
                end
            end
        end 
        function [deltaMz, deltaSz,...
                deltaMzx, deltaSzx] = fcHiddenStateBackwardPassB1(Sz, Sxs,...
                J, mw, deltaM, deltaS, ni, no, gpu) 
            mw  = reshape(mw, [ni, no]);              
            deltaMzx = Sxs;
            deltaSzx = Sxs;
            if isempty(Sxs)
                deltaMloop = repmat(deltaM', [ni, 1]);
                deltaSloop = repmat(deltaS', [ni, 1]);
                if gpu
                    Caz = bsxfun(@times, J, Sz);
                    [deltaMzloop,...
                     deltaSzloop] = arrayfun(@vectorizedDelta_V2, mw, Caz,...
                     deltaMloop, deltaSloop);
                else
                    Caz = J .* Sz;
                    Cwa = mw .* Caz;
                    deltaMzloop = Cwa .* deltaMloop;
                    deltaSzloop = (Cwa.^2) .* deltaSloop;
                end
                deltaMz = sum(deltaMzloop, 2);
                deltaSz = sum(deltaSzloop, 2);
            else                
                deltaMloop = repmat(deltaM', [ni, 1]);
                deltaSloop = repmat(deltaS', [ni, 1]);
                if gpu
                    Caz = bsxfun(@times, J, Sz);
                    Caxs = bsxfun(@times, J, Sxs);
                    [deltaMzloop, deltaSzloop, deltaMxsloop,...
                     deltaSxsloop] = arrayfun(@vectorized4delta, mw, Caz,...
                     Caxs, deltaMloop, deltaSloop);
                else
                    Caz = J .* Sz;
                    Caxs = J .* Sxs;
                    [deltaMzloop, deltaSzloop, deltaMxsloop,...
                     deltaSxsloop] = vectorized4delta(mw, Caz, Caxs,...
                     deltaMloop, deltaSloop);
                end
                deltaMz  = sum(deltaMzloop, 2);
                deltaSz  = sum(deltaSzloop, 2);
                deltaMzx = sum(deltaMxsloop, 2);
                deltaSzx = sum(deltaSxsloop, 2);
            end
        end   
        
        function [mz, Sz] = fcMeanVarMexcuda(mw, Sw, mb, Sb, ma, Sa, widx,...
                bidx, ni, no, B)
            [mz, Sz] = fcMeanVar4matlab(mw, Sw, mb, Sb, ma, Sa, widx, bidx,...
                no, ni, B);
        end
        function [deltaMz, deltaSz,...
                deltaMzx, deltaSzx] = fcHiddenStateBackwardPassMexcuda(Sz,...
                J, mw, deltaM, deltaS, widx, ni, no, B) 
            [deltaMz, deltaSz] = fcHiddenStateBackwardPass4matlab(mw, Sz, J,...
                deltaM, deltaS, widx, ni, no, B);
            deltaMzx = [];
            deltaSzx = [];
        end 
        function [deltaMw, deltaSw,...
                deltaMb, deltaSb] = fcParameterBackwardPassMexcuda(Sw, Sb,...
                ma, deltaM, deltaS, widx, bidx, ni, no, B, deltaMw, deltaSw,...
                deltaMb, deltaSb)  
            [deltaMw, deltaSw, deltaMb,...
             deltaSb] = fcParamBackwardPass4matlab(Sw, ma, Sb, deltaM,...
             deltaS, widx, bidx, ni, B, no, no, B, 1, deltaMw, deltaSw,...
             deltaMb, deltaSb);
        end                      
        
        % Shared functions for update step
        function [deltaM, deltaS] = inovationVector(SzF, dMz, dSz, gpu)
            if gpu
                iSzF  = bsxfun(@rdivide, 1, SzF);
                iSzF(isinf(iSzF)) = zeros(1,1, 'like', dMz);
                [deltaM, deltaS] = arrayfun(@vectorizedDelta, iSzF, dMz, dSz);
            else              
                iSzF   = 1./SzF; 
                iSzF(isinf(iSzF)) = zeros(1,1, 'like', dMz);
                [deltaM, deltaS]  = vectorizedDelta(iSzF, dMz, dSz);
            end           
        end 
        function [deltaMz, deltaSz] = fowardHiddenStateUpdate(mzF, SzF, Cyz,...
                y, gpu)
            if gpu
                dz  = y - mzF;
                SzF = 1 ./ SzF;
                SzF(isinf(SzF)) = 0;
                K = bsxfun(@times, Cyz, SzF);
                deltaMz = bsxfun(@times, K, dz);
                deltaSz = bsxfun(@times, -K, Cyz);
            else
                dz  = y - mzF;
                SzF = 1 ./ SzF;
                SzF(isinf(SzF)) = 0;
                K = Cyz .* SzF;
                deltaMz = K .* dz;
                deltaSz = -K .* Cyz;
            end
        end   
        function theta = globalParameterUpdate(theta, deltaTheta, gpu)          
            [mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx] = tagi.extractParameters(theta);
            [deltaMw, deltaSw, deltaMb, deltaSb, deltaMwx, deltaSwx,...
                deltaMbx, deltaSbx] = tagi.extractParameters(deltaTheta);
            if gpu
                [mw, Sw]   = arrayfun(@twoPlus, mw, Sw, deltaMw, deltaSw);
                [mb, Sb]   = arrayfun(@twoPlus, mb, Sb, deltaMb, deltaSb);
                [mwx, Swx] = arrayfun(@twoPlus, mwx, Swx, deltaMwx, deltaSwx);
                [mbx, Sbx] = arrayfun(@twoPlus, mbx, Sbx, deltaMbx, deltaSbx);
            else
                [mw, Sw]   = twoPlus(mw, Sw, deltaMw, deltaSw);
                [mb, Sb]   = twoPlus(mb, Sb, deltaMb, deltaSb);
                [mwx, Swx] = twoPlus(mwx, Swx, deltaMwx, deltaSwx);
                [mbx, Sbx] = twoPlus(mbx, Sbx, deltaMbx, deltaSbx);
            end
            theta = tagi.compressParameters(mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx);
        end
        function theta = globalParameterUpdateMultiGPUs(theta, deltaTheta,...
                numParamsPerlayer, numDevices)  
            numParams  = sum(numParamsPerlayer, 2);           
            deltaTheta = cat(2, deltaTheta{:});
            [deltaMw, deltaSw, deltaMb, deltaSb, deltaMwx, deltaSwx,...
                deltaMbx, deltaSbx] = tagi.extractParameters_V2(deltaTheta);
            deltaMw  = cat(1, deltaMw{:});
            deltaSw  = cat(1, deltaSw{:});
            deltaMb  = cat(1, deltaMb{:});
            deltaSb  = cat(1, deltaSb{:});
            deltaMwx = cat(1, deltaMwx{:});
            deltaSwx = cat(1, deltaSwx{:});
            deltaMbx = cat(1, deltaMbx{:});
            deltaSbx = cat(1, deltaSbx{:});  
            
            deltaMw  = sum(reshape(cat(1, deltaMw{:}), [numParams(1), numDevices]), 2);
            deltaSw  = sum(reshape(cat(1, deltaSw{:}), [numParams(1), numDevices]), 2);
            deltaMb  = sum(reshape(cat(1, deltaMb{:}), [numParams(2), numDevices]), 2);
            deltaSb  = sum(reshape(cat(1, deltaSb{:}), [numParams(2), numDevices]), 2);
            deltaMwx = sum(reshape(cat(1, deltaMwx{:}), [numParams(3), numDevices]), 2);
            deltaSwx = sum(reshape(cat(1, deltaSwx{:}), [numParams(3), numDevices]), 2);
            deltaMbx = sum(reshape(cat(1, deltaMbx{:}), [numParams(4), numDevices]), 2);
            deltaSbx = sum(reshape(cat(1, deltaSbx{:}), [numParams(4), numDevices]), 2);            
            [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = tagi.extractParameters(theta);
            [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = tagi.catParameters(mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx);            
            [mw, Sw]   = arrayfun(@twoPlus, mw, Sw, deltaMw, deltaSw);
            [mb, Sb]   = arrayfun(@twoPlus, mb, Sb, deltaMb, deltaSb);
            [mwx, Swx] = arrayfun(@twoPlus, mwx, Swx, deltaMwx, deltaSwx);
            [mbx, Sbx] = arrayfun(@twoPlus, mbx, Sbx, deltaMbx, deltaSbx);
            [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = tagi.distributeParameters2Layers(mw, Sw, mb, Sb,...
                mwx, Swx, mbx, Sbx, numParamsPerlayer);
            theta = tagi.compressParameters(mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx);
        end
        function [mnorm, Snorm] = distrNorm(m, S)
            N     = length(S);
            mhat  = mean(m);
            Shat  = 1 / N * (sum(S, 1)+ sum((m - mhat).^2, 1)) + 1E-8;
            mnorm = (m - mhat) ./ (sqrt(abs(Shat)));
            Snorm = (1 ./ abs(Shat)) .* S;
        end         
        function theta = globalParameterUpdateCUDA(theta, deltaTheta,...
                wxupdate, numParams)          
            [mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx] = tagi.extractParameters(theta);
            [deltaMw, deltaSw, deltaMb, deltaSb, deltaMwx, deltaSwx,...
                deltaMbx, deltaSbx] = tagi.extractParameters(deltaTheta);
            
            [mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx] = globalParamUpdate(mw, Sw,...
                mb, Sb, mwx, Swx, mbx, Sbx, deltaMw, deltaSw, deltaMb, deltaSb,...
                deltaMwx, deltaSwx, deltaMbx, deltaSbx, numParams(1, end),...
                numParams(2, end), numParams(3, end), numParams(4, end), wxupdate);
            
            theta = tagi.compressParameters(mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx);
        end
        
        % Noise update
        function [l, v2] = detachMeanVar(x, nl, nv2, B, rB)
            x  = reshape(x, [nl + nv2, B * rB]);
            l  = reshape(x(1 : nl, :), [B * nl, rB]);
            v2 = reshape(x(nl + 1 : end, :), [B * nv2, rB]);
        end
        function x = attachMeanVar(l, v2, nl, nv2, B, rB)
            l  = reshape(l, [nl, B * rB]);
            v2 = reshape(v2, [nv2, B * rB]);
            x  = [l; v2];
            x  = reshape(x, [(nl + nv2) * B, rB]);
        end
        function [deltaMlz, deltaSlz, deltaMv2z,...
                deltaSv2z] = noiseUpdate4regression(Slz, mla, Sla, J, Jv2,...
                mv2a, Sv2a, Cv2a, y, Sy, sv, net)
            Cyz = J .* Slz ;
            Syf = Sla + mv2a + sv.^2;
            
            if net.imperfect_obs
                deltaMz = y - mla;
                deltaSz = Sy - Syf;
                % Updating the hidden node(s) used for estimating y
                [deltaMlz, deltaSlz] = tagi.inovationVector(Syf,...
                    deltaMz, deltaSz, net.gpu);
                % RTS update
                deltaMlz = Sla .* deltaMlz;
                deltaSlz = Sla .* deltaSlz .* Sla;
                % Updating the hidden node(s) used for estimating the noise
                [deltaMv, deltaSv] = tagi.inovationVector(mv2a,...
                    deltaMz, deltaSz, net.gpu);
                % RTS update
                deltaMv = mv2a .* deltaMv;
                deltaSv = mv2a .* deltaSv .* mv2a;
            else
                [deltaMlz, deltaSlz] = tagi.fowardHiddenStateUpdate(mla, Syf,...
                    Cyz, y, net.gpu);
                [deltaMv, deltaSv] = tagi.fowardHiddenStateUpdate(mla, Syf,...
                    mv2a, y, net.gpu);
            end

            mvUd = deltaMv;
            SvUd = mv2a + deltaSv;
            
            % Update activated standard deviation for z
            yv2  = mvUd.^2 + SvUd;
            Sv2f = 2 * SvUd.^2 + 4 * (mvUd.^2) .* SvUd;   
            [deltaMv2z, deltaSv2z] = tagi.noiseBackwardUpdate(mv2a,...
                3 * Sv2a + 2 * mv2a.^2, Jv2 .* Cv2a, yv2, Sv2f, net.gpu); 
        end
        function [deltaMlz, deltaSlz, deltaMv2z,...
                deltaSv2z] = homoNoiseUpdate4regression(Slz, mla, Sla, J,...
                mv2a, Sv2a, y, gpu)
            Cyz = J .* Slz ;
            Syf = Sla + mv2a;
            
            [deltaMlz, deltaSlz] = tagi.fowardHiddenStateUpdate(mla, Syf,...
                Cyz, y, gpu);
            [deltaMv, deltaSv] = tagi.fowardHiddenStateUpdate(mla, Syf,...
                mv2a, y, gpu);
            
            mvUd = deltaMv;
            SvUd = mv2a + deltaSv;
            
            % Update activated standard deviation for z
            yv2  = mvUd.^2 + SvUd;
            Sv2f = 2 * SvUd.^2 + 4 * (mvUd.^2) .* SvUd; 
            Cv2a = Sv2a;
            [deltaMv2z, deltaSv2z] = tagi.noiseBackwardUpdate(mv2a,...
                3 * Sv2a + 2 * mv2a.^2, Cv2a, yv2, Sv2f, gpu); 
        end
        function [deltaMz, deltaSz] = noiseBackwardUpdate(maF, SaF, CzzF,...
                maB, SaB, gpu)            
            if gpu
                funM    = @(x, y, z) x .* (y - z);
                funS    = @(x, y, z) x .* (y - z) .* x;
                Jz      = CzzF ./ SaF; 
                deltaMz = arrayfun(funM, Jz, maB, maF);
                deltaSz = arrayfun(funS, Jz, SaB, SaF);
            else
                Jz      = CzzF ./ SaF; 
                deltaMz = Jz .* (maB - maF);
                deltaSz = Jz .* (SaB - SaF) .* Jz;
            end
        end  
        
        % Initialization for weights and bias   
        function theta = initializeWeightBias(net)
            %Initialization
            nodes     = double(net.nodes);
            numLayers = length(net.nodes);
            layer     = net.layer;
            idxw      = net.idxw;
            idxwXsc   = net.idxwXsc;
            idxbXsc   = net.idxbXsc;
            idxb      = net.idxb;
            biasStd   = 1E-2;
            gainMw    = cast(net.gainMw, net.dtype);
            gainSw    = cast(net.gainSw, net.dtype); 
            gainMb    = cast(net.gainMb, net.dtype);
            gainSb    = cast(net.gainSb, net.dtype); 
            mw        = tagi.createInitCellwithArray(numLayers-1, net.dtype,...
                net.gpu);
            Sw        = mw;
            mb        = mw;
            Sb        = mw;
            mwx       = mw;
            Swx       = mw;
            mbx       = mw;
            Sbx       = mw;
            if net.gain_v2
                gainMv2 = cast(net.gainMv2, net.dtype);
                gainSv2 = cast(net.gainSv2, net.dtype);
            end
            for j = 2:numLayers
                if ~isempty(idxw{j-1})                    
                    if layer(j) == net.layerEncoder.conv ...
                            || layer(j) == net.layerEncoder.tconv % Conv. layer
                        fanIn  = (cast(net.kernelSize(j-1), net.dtype) .^2) ...
                            * cast(net.filter(j-1), net.dtype);
                        if net.xsc(j-1)~=0
                            fanIn = 2 * fanIn;
                        end
                        if strcmp(net.initParamType, 'Xavier')
                            if j < numLayers ...
                                    && (layer(j+1) == net.layerEncoder.mp ...
                                    || layer(j+1) == net.layerEncoder.ap)
                                fanOut = ((cast(net.kernelSize(j-1), net.dtype) .^ 2) ...
                                    * cast(net.filter(j), net.dtype)) ...
                                    / (cast(net.kernelSize(j), net.dtype) .^2);
                            else
                                fanOut = ((cast(net.kernelSize(j-1), net.dtype).^2) ...
                                    * cast(net.filter(j), net.dtype));
                            end
                            scale = 2 / (fanIn + fanOut);
                            Sw{j-1} = (gainSw(j-1)) * scale * ...
                                ones(length(idxw{j-1}), 1, net.dtype);
                        elseif strcmp(net.initParamType, 'He')
                            scale   = 1 / (fanIn);
                            Sw{j-1} = (gainSw(j-1)) * scale * ...
                                ones(length(idxw{j-1}), 1, net.dtype);
                        end 
                        mw{j-1} = gainMw(j-1) * randn(length(Sw{j-1}), 1) ...
                            .* sqrt(Sw{j-1});
                        if ~isempty(idxb{j-1})
                            Sb{j-1} = gainSb(j-1) * biasStd * ...
                                ones(length(idxb{j-1}), 1, net.dtype);
                            mb{j-1} = gainMb(j-1) * randn(length(Sb{j-1}), 1) ...
                                .* sqrt(Sb{j-1});
                        end
                    elseif layer(j) == net.layerEncoder.ln ...
                            || layer(j) == net.layerEncoder.bn
                        Sb{j-1} = 1E-4 * gainSw(j-1) * ...
                            ones(length(idxb{j-1}), 1, net.dtype);
                        mb{j-1} = 0 * rand(length(Sb{j-1}), 1, net.dtype) ...
                            .* sqrt(Sb{j-1});
                        Sw{j-1} = 1 * ones(length(idxw{j-1}), 1, net.dtype);
                        mw{j-1} = 1 * ones(length(idxw{j-1}), 1, net.dtype);
                    else
                        fanIn  = nodes(j-1);
                        fanOut = nodes(j);
                        if strcmp(net.initParamType, 'Xavier')
                            scale = 2/(fanIn + fanOut);
                            Sw{j-1} = (gainSw(j-1)) * scale * ...
                                ones(length(idxw{j-1}), 1, net.dtype);
                        elseif strcmp(net.initParamType, 'He')
                            scale = 1 / fanIn;
                            Sw{j-1} = (gainSw(j-1)) * scale * ...
                                ones(length(idxw{j-1}), 1, net.dtype);
                            if net.gain_v2 && j == numLayers
                                Sw{j-1}(fanIn*(fanOut-net.nLchol)+1:end) = gainSv2 * scale * ...
                                    ones(length(idxw{j-1}(fanIn*(fanOut-net.nLchol)+1:end)), 1, net.dtype); 
                            end
                        end
                        mw{j-1} = gainMw(j-1)*randn(length(Sw{j-1}), 1) ...
                            .* sqrt(Sw{j-1});
                        if net.gain_v2
                        mw{j-1}(fanIn*(fanOut-net.nLchol)+1:end) = gainMv2*randn(length(Sw{j-1}(fanIn*(fanOut-net.nLchol)+1:end)), 1) ...
                            .* sqrt(Sw{j-1}(fanIn*(fanOut-net.nLchol)+1:end));
                        end
                        if ~isempty(idxb{j-1})
                            Sb{j-1} = gainSb(j-1) * scale * ...
                                ones(length(idxb{j-1}), 1, net.dtype);
                            mb{j-1} = gainSb(j-1) * randn(length(Sb{j-1}), 1) ...
                                .* sqrt(Sb{j-1});
                        end
                    end  
                end 
                if net.xsc(j) ~=0 ...
                        && (net.filter(net.xsc(j)) ~= net.filter(j) ...
                        || net.imgW(net.xsc(j)) ~= net.imgW(j))
                    idxXsc = net.xsc(j);                                     
                    fanIn  = cast(net.filter(idxXsc), net.dtype);
                    fanOut = cast(net.filter(j), net.dtype);
                    if strcmp(net.initParamType, 'Xavier')
                        Swx{idxXsc} = (gainS(idxXsc)) * (2 / (fanIn + fanOut)) ...
                            * ones(length(idxwXsc{idxXsc}), 1, net.dtype);
                    elseif strcmp(net.initParamType, 'He')
                        Swx{idxXsc} = (1 / (fanIn)) * ...
                            ones(length(idxwXsc{idxXsc}), 1, net.dtype);
                    end
                    mwx{idxXsc} = randn(length(Swx{idxXsc}), 1) ...
                        .* sqrt(Swx{idxXsc});
                    if ~isempty(idxbXsc{idxXsc})
                        Sbx{idxXsc} = 1E-6 * ones(length(idxbXsc{idxXsc}), 1, net.dtype);
                        mbx{idxXsc} = 0 * randn(length(Sbx{idxXsc}), 1) ...
                            .* sqrt(Sbx{idxXsc});
                    end                   
                    if net.gpu
                        mwx{idxXsc} = gpuArray(mwx{idxXsc});
                        Swx{idxXsc} = gpuArray(Swx{idxXsc});
                        mbx{idxXsc} = gpuArray(mbx{idxXsc});
                        Sbx{idxXsc} = gpuArray(Sbx{idxXsc});
                    end
                end
                clear fanIn
                % Send to gpu
                if net.gpu
                    mw{j-1} = gpuArray(mw{j-1});
                    Sw{j-1} = gpuArray(Sw{j-1});
                    mb{j-1} = gpuArray(mb{j-1});
                    Sb{j-1} = gpuArray(Sb{j-1});                    
                end
            end 
            [mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx] = tagi.catParameters(mw, Sw,...
                mb, Sb, mwx, Swx, mbx, Sbx);
           theta = tagi.compressParameters(mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx); 
        end      
        function states = initializeStates(nodes, B, rB, xsc, dtype, gpu)
            % Normal net
            numLayers = length(nodes);          
            mz  = tagi.createStateCellarray(nodes, numLayers, B, rB, 0, dtype, gpu); 
            Sz  = mz; 
            ma  = mz;
            Sa  = mz;
            J   =  tagi.createStateCellarray(nodes, numLayers, B, rB, 1, dtype, gpu);
            % Residual net
            idx = xsc~=0;
            mdxs = cell(numLayers, 1);
            mdxs(idx) = mz(idx);
            Sdxs = mdxs;
            mxs  = mdxs;
            Sxs  = mdxs;
            states = tagi.compressStates(mz, Sz, ma, Sa, J, mdxs, Sdxs, mxs, Sxs);
        end
        function [deltaxs, deltadxs] = initializeShortcutStateDelta(xsc,...
                idxXsc, x, B, rB)
            layers   = xsc(xsc~=0);
            deltaxs  = cell(length(xsc), 1);
            deltadxs = cell(length(xsc), 1);
            for j = layers
                if ~isempty(idxXsc{j})
                    deltaxs{j}  = zeros(length(idxXsc{j})*B, rB, 'like', x{j});
                    deltadxs{j} = deltaxs{j};
                else
                    deltadxs{j} = zeros(size(x{j}), 'like', x{j});
                    deltaxs{j}  = zeros(size(x{j}), 'like', x{j});
                end
            end
        end
        function states = initializeInputs(states, mz0, Sz0, ma0, Sa0, J0,...
                mdxs0, Sdxs0, mxs0, Sxs0, xsc)
            [mz, Sz, ma, Sa, J, mdxs, Sdxs, mxs, Sxs] = tagi.extractStates(states);
            % Normal net
            mz{1} = mz0;
            if any(isempty(Sz0))
                Sz{1} = zeros(size(mz0), 'like', mz0);
            else
                Sz{1} = Sz0;
            end
            if any(isempty(ma0))
                ma{1} = mz0;
            else
                ma{1} = ma0;
            end 
            if any(isempty(Sa0))
                Sa{1} = Sz{1};
            else
                Sa{1} = Sa0;
            end   
            if any(isempty(J0))
                J{1} = ones(size(mz0), 'like', mz0);
            else
                J{1} = J0;
            end  
            % Residual net
            if any(isempty(mdxs0))&&~all(xsc==0)
                mdxs{1} = mz0;
            else
                mdxs{1} = mdxs0;
            end
            if any(isempty(Sdxs0))&&~all(xsc==0)
                Sdxs{1} = zeros(size(mz0), 'like', mz0);
            else
                Sdxs{1} = Sdxs0;
            end
            if any(isempty(mxs0))&&~all(xsc==0)
                mxs{1} = mz0;
            else
                mxs{1} = mxs0;
            end
            if any(isempty(Sxs0))&&~all(xsc==0)
                Sxs{1} = zeros(size(mz0), 'like', mz0);
            else
                Sxs{1} = Sxs0;
            end
            states = tagi.compressStates(mz, Sz, ma, Sa, J, mdxs, Sdxs, mxs, Sxs);
        end
        function maxIdx = initializeMaxPoolingIndices(nodes, layers,...
                layerEncoder, B, rB, dtype, gpu)
            if gpu
                zeroPad = zeros(1, 1, dtype, 'gpuArray');
            else
                zeroPad = zeros(1, 1, dtype);
            end
            numLayers = length(nodes);
            maxIdx = cell(numLayers, 1);
            maxPoolingLayers = find(layers==layerEncoder.mp);
            if ~isempty(maxPoolingLayers)
                for j = maxPoolingLayers
                    maxIdx{j} = zeros(nodes(j)*B, rB, 'like', zeroPad);
                end
            end
        end
        function normStat = initializeNormStat(nodes, filter, B, rB, layers,...
                layerEncoder, x)
            numLayers = length(nodes);
            mra = cell(numLayers, 1);
            layNorm = layers==layerEncoder.ln;
            batNormConv = layers==layerEncoder.bn&(layers==layerEncoder.conv|layers==layerEncoder.tconv|layers==layerEncoder.mp|layers==layerEncoder.ap);
            batNormfc = layers==layerEncoder.bn&layers==layerEncoder.fc;
            for j = layNorm
                mra{j} = zeros(B, rB, 'like', x);
            end
            for j = batNormfc
                mra{j} = zeros(nodes(j), rB, 'like', x);
            end
            for j = batNormConv
                mra{j} = zeros(filter(j), rB, 'like', x);
            end
            Sra = mra;
            normStat = tagi.compressNormStat(mra, Sra);
        end  
        function deltaTheta = initializeDeltaTheta(theta, rB, numLayers)
            deltaTheta = cell(numLayers-1, 1);
            for j = 1:numLayers-1
                deltaTheta{j} = repmat(theta{j}, [1, rB]);
            end
        end
        function [mw, Sw, mb, Sb, mwx, Swx, mbx, Sbx] = catParameters(mw,...
                Sw, mb, Sb, mwx, Swx, mbx, Sbx)
            mw  = cat(1, mw{:});
            Sw  = cat(1, Sw{:});
            mb  = cat(1, mb{:});
            Sb  = cat(1, Sb{:});
            mwx = cat(1, mwx{:});
            Swx = cat(1, Swx{:});
            mbx = cat(1, mbx{:});
            Sbx = cat(1, Sbx{:});
        end
        function [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = distributeParameters2Layers(mw, Sw, mb, Sb, mwx,...
                Swx, mbx, Sbx, numParams)
            mw  = mat2cell(mw, numParams(1, :));
            Sw  = mat2cell(Sw, numParams(1, :));
            mb  = mat2cell(mb, numParams(2, :));
            Sb  = mat2cell(Sb, numParams(2, :));
            mwx = mat2cell(mwx, numParams(3, :));
            Swx = mat2cell(Swx, numParams(3, :));
            mbx = mat2cell(mbx, numParams(4, :));
            Sbx = mat2cell(Sbx, numParams(4, :));
        end           
               
        % Storing
        function states = compressStates(mz, Sz, ma, Sa, J, mdxs, Sdxs, mxs,...
                Sxs)
            states = cell(9, 1);
            states{1} = mz;
            states{2} = Sz;
            states{3} = ma;
            states{4} = Sa;
            states{5} = J;
            states{6} = mdxs;
            states{7} = Sdxs;
            states{8} = mxs;
            states{9} = Sxs;
        end
        function [mz, Sz, ma, Sa, J, mdxs, Sdxs,...
                mxs, Sxs] = extractStates(states)
            mz   = states{1};
            Sz   = states{2};
            ma   = states{3};
            Sa   = states{4};
            J    = states{5};
            mdxs = states{6};
            Sdxs = states{7};
            mxs  = states{8};
            Sxs  = states{9};
        end
        function [mz, Sz, ma, Sa, J, mdxs, Sdxs,...
                mxs, Sxs] = extractStatesMultiGPUs(states)
            spmd
                mz   = states{1};
                Sz   = states{2};
                ma   = states{3};
                Sa   = states{4};
                J    = states{5};
                mdxs = states{6};
                Sdxs = states{7};
                mxs  = states{8};
                Sxs  = states{9};
            end
        end
        function theta = compressParameters(mw, Sw, mb, Sb, mwx, Swx, mbx,...
                Sbx)
            theta     = cell(8, 1);
            theta{1}  = mw;
            theta{2}  = Sw;
            theta{3}  = mb;
            theta{4}  = Sb;
            theta{5}  = mwx;
            theta{6}  = Swx;
            theta{7}  = mbx;
            theta{8}  = Sbx;
        end
        function [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = extractParameters(theta)
            mw  = theta{1};
            Sw  = theta{2};
            mb  = theta{3};
            Sb  = theta{4};
            mwx = theta{5};
            Swx = theta{6};
            mbx = theta{7};
            Sbx = theta{8};
        end
        function [mw, Sw, mb, Sb, mwx, Swx,...
                mbx, Sbx] = extractParameters_V2(theta)
            mw  = theta(1, :);
            Sw  = theta(2, :);
            mb  = theta(3, :);
            Sb  = theta(4, :);
            mwx = theta(5, :);
            Swx = theta(6, :);
            mbx = theta(7, :);
            Sbx = theta(8, :);
        end
        function normStat = compressNormStat(mra, Sra)
            normStat = cell(2, 1);
            normStat{1} = mra;
            normStat{2} = Sra;
        end
        function [mra, Sra] = extractNormStat(normStat)
            mra = normStat{1};
            Sra = normStat{2};
        end   
        
        % Create cell with an array
        function x = createInitCellwithArray(numLayers, dtype, gpu)
            x = cell(numLayers, 1);
            if gpu
                x(:) = {gpuArray(nan(1, 1, dtype))};
            else
                x(:) = {nan(1, 1, dtype)};
            end
        end
        function z = createStateCellarray(nodes, numLayers, B, rB, value,...
                dtype, gpu)   
            z = cell(numLayers, 1);
            if gpu
                zeroPad = zeros(1,1,dtype, 'gpuArray');
            else
                zeroPad = zeros(1,1,dtype);
            end
            for j = 2:numLayers               
                z{j} = zeros(nodes(j)*B, rB, 'like', zeroPad) + value;
            end
        end 
        function d = createDevCellarray(nodes, numLayers, B, rB, dtype, gpu)   
            d = cell(numLayers, 1);
            if gpu
                onePad = ones(1,1,dtype, 'gpuArray');
            else
                onePad = ones(1,1,dtype);
            end
            for j = 1:numLayers               
                d{j} = ones(nodes(j)*B, rB, 'like', onePad);
            end
%             d{numLayers} = zeroPad + 1;
        end
        function normStat = createInitNormStat(net)
            mra   = cell(length(net.nodes) -1, 1);            
            Sra   = cell(length(net.nodes) -1, 1);
            dtype = net.dtype;
            numLayers = length(net.layer);
            layer     = net.layer;
            if net.gpu
                zeroPad = zeros(1,1,dtype, 'gpuArray');
            else
                zeroPad = zeros(1,1,dtype);
            end
            for j = 2:numLayers 
                if layer(j-1) == net.layerEncoder.conv
                    if layer(j) == net.layerEncoder.ln
                        mra{j-1} = zeros(net.batchSize * net.repBatchSize, 1, 'like', zeroPad);
                        Sra{j-1} = ones(net.batchSize * net.repBatchSize, 1, 'like', zeroPad);
                    elseif layer(j) == net.layerEncoder.bn
                        mra{j-1} = zeros(net.filter(j), 1, 'like', zeroPad);
                        Sra{j-1} = ones(net.filter(j), 1, 'like', zeroPad);
                    else
                        mra{j-1} = zeros(1, 1, 'like', zeroPad);
                        Sra{j-1} = ones(1, 1, 'like', zeroPad);
                    end
                elseif layer(j-1) == net.layerEncoder.fc
                    if layer(j) == net.layerEncoder.ln
                        mra{j-1} = zeros(net.batchSize * net.repBatchSize, 1, 'like', zeroPad);
                        Sra{j-1} = ones(net.batchSize * net.repBatchSize, 1, 'like', zeroPad);
                    elseif layer(j) == net.layerEncoder.bn
                        mra{j-1} = zeros(net.nodes(j), 1, 'like', zeroPad);
                        Sra{j-1} = ones(net.nodes(j), 1, 'like', zeroPad);
                    else
                        mra{j-1} = zeros(1, 1, 'like', zeroPad);
                        Sra{j-1} = ones(1, 1, 'like', zeroPad);
                    end
                end
            end
            normStat = tagi.compressNormStat(mra, Sra);
        end
    end
end