%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         indices
% Description:  Build indices for Tractable Approximate Gaussian Inference
% Authors:      Luong-Ha Nguyen & James-A. Goulet
% Created:      November 3, 2019
% Updated:      August 27, 2021
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef indices
    methods (Static)
        function netInfos = savedInfo(net)
            if strcmp(net.task, 'classification')...
                    || strcmp(net.task, 'discrimination') ...
                    || strcmp(net.task, 'generation')...
                    || strcmp(net.task, 'autoencoding')...
                    || strcmp(net.task, 'encoding')...
                    || strcmp(net.task, 'decoding')...
                    ||strcmp(net.task, 'sharing')
                
                netInfos.filter      = net.filter;
                netInfos.kernelSize  = net.kernelSize;
                netInfos.padding     = net.padding;
                netInfos.paddingType = net.paddingType;
                netInfos.stride      = net.stride;
                netInfos.actFunIdx   = net.actFunIdx;
                netInfos.imgSize     = net.imgSize;
                netInfos.svDecayFactor = net.svDecayFactor;
                netInfos.svmin       = net.svmin;
                netInfos.da          = net.da;
                netInfos.imgW        = net.imgW;
                netInfos.imgH        = net.imgH;
            end
            netInfos.layer         = net.layer;
            netInfos.actFunIdx     = net.actFunIdx; 
            netInfos.nodes         = net.nodes;
            netInfos.sv            = net.sv;
            netInfos.batchSize     = net.batchSize;
            netInfos.repBatchSize  = net.repBatchSize;
            netInfos.initParamType = net.initParamType;
            netInfos.gainMw        = net.gainMw;
            netInfos.gainSw        = net.gainSw;
            netInfos.gainMb        = net.gainMb;
            netInfos.gainSb        = net.gainSb;
            netInfos.seed          = net.seed;
            netInfos.svmin         = net.svmin;
            if isfield(net, 'rl')
                netInfos.rl = net.rl;
            end
        end
        function net = initialization(net)
            numLayers = length(net.layer);
            if strcmp(net.task, 'classification')...
                    || strcmp(net.task, 'discrimination') ...
                    || strcmp(net.task, 'generation')...
                    || strcmp(net.task, 'autoencoding')...
                    || strcmp(net.task, 'encoding')...
                    || strcmp(net.task, 'decoding')...
                    ||strcmp(net.task, 'sharing')...
                    ||strcmp(net.task, 'regression')
                if ~isfield(net, 'filter')
                    net.filter = zeros(1, numLayers);
                end
                if ~isfield(net, 'kernelSize')
                    net.kernelSize = zeros(1, numLayers);
                end
                if ~isfield(net, 'padding')
                    net.padding = zeros(1, numLayers);
                end
                if ~isfield(net, 'paddingType')
                    net.paddingType = zeros(1, numLayers);
                end
                if ~isfield(net, 'stride')
                    net.stride = zeros(1, numLayers);
                end
                if ~isfield(net, 'actFunIdx')
                    net.actFunIdx = zeros(1, numLayers);
                end
                if ~isfield(net, 'actBound')
                    net.actBound = ones(1, numLayers);
                end                                              
                if ~isfield(net, 'imgH')
                    net.imgH = zeros(1, numLayers);
                    if isfield(net, 'imgSize')
                        net.imgH(1) = net.imgSize(1);
                    end
                end
                if ~isfield(net, 'imgW')
                    net.imgW = zeros(1, numLayers);
                    if isfield(net, 'imgSize')
                        net.imgW(1) = net.imgSize(1);
                    end
                end
                if ~isfield(net, 'obsShow')
                    net.obsShow = 10000;
                end
                if ~isfield(net, 'epsilon')
                    net.epsilon = 1E-4;
                end
                if ~isfield(net, 'normMomentum')
                    net.normMomentum = 0.0;
                    net.normMomentumRef = 0.9;
                end
                if ~isfield(net, 'da')
%                     net = indices.daTypeEncoder(net);
                    net.da.enable = 0;
                    net.da.p      = 0.5;
                    net.da.types  = [];
                end                
                if ~isfield(net, 'earlyStop')
                    net.earlyStop = 0;
                end
                if ~isfield(net, 'displayMode')
                    net.displayMode = 1;
                end 
                if ~isfield(net, 'errorRateEval')
                    net.errorRateEval = 1;
                end
                if ~isfield(net, 'numDevices')
                    net.numDevices = 1;
                end
            end
            
            % Resnet
            if ~isfield(net, 'xsc')
                net.xsc = zeros(1, numLayers);
                % first defaut short cut (see hiddenStateBackpass)
                %                     net.xsc(1) = 1;
            end
            if ~isfield(net, 'wxupdate')
                net.wxupdate = false;
            end
            
            % Parameters
            if ~isfield(net, 'initParamType')
                net.initParamType = 'Xavier';
            end
            if ~isfield(net, 'gainMw')
                net.gainMw = 1*ones(1, numLayers-1);
            end
            if ~isfield(net, 'gainSw')
                if strcmp(net.initParamType, 'Xavier')
                    net.gainSw = 1*ones(1, numLayers-1);
                elseif strcmp(net.initParamType, 'He')
                    net.gainSw = 1*ones(1, numLayers-1);
                end
            end
            if ~isfield(net,'gainMb')
                net.gainMb = 1*ones(1, numLayers-1);
            end
            if ~isfield(net, 'gainSb')
                if strcmp(net.initParamType, 'Xavier')
                    net.gainSb = 1*ones(1, numLayers-1);
                elseif strcmp(net.initParamType, 'He')
                    net.gainSb = 1*ones(1, numLayers-1);
                end
            end
            if ~isfield(net, 'maxEpoch')
                net.maxEpoch = 10;
            end                       
            if ~isfield(net, 'imgSize')
                net.imgSize = [0 0 0];
            end
            if ~isfield(net, 'savedEpoch')
                net.savedEpoch = 0;
            end                        
            if ~isfield(net,'seed')
                net.seed = nan;
            end
            
            % Activation
            if ~isfield(net, 'lreluRate')
                net.lreluRate = 0.2;
            end
            
            % Input
            if ~isfield(net, 'convariateEstm')
                net.convariateEstm = false;
            end
            if ~isfield(net, 'fullcov')
                net.fullcov = false;
            end
            if ~isfield(net, 'PO')
                net.PO = false;
            end            
            
            % Outputs
            if ~isfield(net,'lastLayerUpdate')
                net.lastLayerUpdate = true;
            end
            if ~isfield(net, 'learnSv')
                net.learnSv = 0;
            end
            if ~isfield(net, 'nye')
                net.nye = net.ny;
            end
            if ~isfield(net, 'isUdidx')
                net.isUdidx = false;
            end
            if ~isfield(net, 'svDecayFactor')
                net.svDecayFactor = 0.8;
            end
            if ~isfield(net, 'svmin')
                net.svmin = 0.0;
            end
            if ~isfield(net, 'learningRateSchedule')
                net.learningRateSchedule = 0;
            end
            if ~isfield(net, 'scheduledSv')
                net.scheduledSv = 0;
            end
            if ~isfield(net, 'imperfect_obs')
                net.imperfect_obs = false;
            end
            
            % Derivative
            if ~isfield(net, 'collectDev')
                net.collectDev = false;
            end
            
            % C/CUDA compiler & GPU
            if ~isfield(net, 'cuda')
                net.cuda = false;
            end
            if ~isfield(net, 'gpu')
                net.gpu = false;
            end
                      
        end

        function net = layerEncoder(net)
            % Full-conetected layer
            layerEncoder.fc   = 1*ones(1, 1);
            % Convolutional layer
            layerEncoder.conv = 2*ones(1, 1);
            % Transposed convolutational layer
            layerEncoder.tconv = 21*ones(1, 1);
            % Max pooling layer 
            layerEncoder.mp   = 3*ones(1, 1);
            layerEncoder.mup  = 31*ones(1, 1);
            % Average pooling layer 
            layerEncoder.ap   = 4*ones(1, 1); 
            % Layer normalization
            layerEncoder.ln   = 5*ones(1, 1);
            % Batch Normalization
            layerEncoder.bn   = 6*ones(1, 1);                      
            % Output
            net.layerEncoder   = layerEncoder;           
        end
        function net = parameters(net)
            % See document for the parameter's ordering
            % Initialization   
            
            if strcmp(net.dtype, 'single')
                dtype = 'single';
            elseif strcmp(net.dtype, 'double')
                dtype = 'double';
            end
            if net.cuda; dtype = 'double'; end
            
            net.nodes     = cast(net.nodes, dtype);
            net.batchSize = cast(net.batchSize, dtype);
            net.imgSize   = cast(net.imgSize, dtype);
            net.kernelSize= cast(net.kernelSize, dtype);
            net.filter    = cast(net.filter, dtype);
            net.stride    = cast(net.stride, dtype);
            net.imgW      = cast(net.imgW, dtype);
            net.imgH      = cast(net.imgH, dtype);
            
            nodes     = net.nodes;
            layer     = net.layer;
            numLayers = length(nodes);           
            % Bias
            idxb      = cell(numLayers - 1, 1);
            % Weights
            idxw      = cell(numLayers - 1, 1);
            idxwXsc   = cell(numLayers - 1, 1);
            idxbXsc   = cell(numLayers - 1, 1);
            
            % Total number of parameters
            numParams         = cell(numLayers - 1, 1);
            numParamsPerlayer = ones(4, numLayers - 1, 'like', nodes);
            paramUpdateIdx    = zeros(4, numLayers - 1);
            totalNumParams    = 0;
            for j = 1:numLayers-1
                if j > 1 && j < numLayers%-sum(layer==1)+1
                    if layer(j) == net.layerEncoder.conv ...
                            || layer(j) == net.layerEncoder.tconv ...
                            || layer(j) == net.layerEncoder.mp ...
                            || layer(j) == net.layerEncoder.ap
                        nodes(j) = net.imgW(j) * net.imgH(j) * net.filter(j);
                    elseif layer(j) == net.layerEncoder.ln ...
                            || layer(j) == net.layerEncoder.bn
                        nodes(j) = nodes(j-1);
                    end
                end
                if layer(j+1) == 1 % Full conetected layer
                    if j < numLayers-2 ...
                            && (layer(j+2) == net.layerEncoder.ln ...
                            || layer(j+2) == net.layerEncoder.bn)
                        numParams{j} = nodes(j+1)*nodes(j);
                    else
                        numParams{j} = nodes(j+1) * nodes(j) + nodes(j+1);
                        idxb{j} = colon(1, nodes(j+1))';                      
                    end
                    idxw{j} = colon(1, nodes(j+1) * nodes(j))';
                elseif layer(j+1) == net.layerEncoder.conv % Conv layer
                    if net.paddingType(j) == 1
                        imgWloop = (double(net.imgW(j)) -...
                            double(net.kernelSize(j)) +...
                            2 * double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                        imgHloop = (double(net.imgH(j)) -...
                            double(net.kernelSize(j)) +...
                            2 * double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                    else
                        imgWloop = (double(net.imgW(j)) -...
                            double(net.kernelSize(j)) +...
                            double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                        imgHloop = (double(net.imgH(j)) -...
                            double(net.kernelSize(j)) +...
                            double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                    end
                    if floor(imgWloop) ~= imgWloop ...
                            || floor(imgHloop) ~= imgHloop
                        error('The hyperparameters for conv. layer are invalid')
                    else
                        net.imgW(j+1) = imgWloop;
                        net.imgH(j+1) = imgHloop;
                    end
                    if j < numLayers-2 ...
                            && (layer(j+2) == net.layerEncoder.ln ...
                            || layer(j+2) == net.layerEncoder.bn) 
                        numParams{j} = net.kernelSize(j) * net.kernelSize(j)...
                            * net.filter(j) * net.filter(j+1);
%                         numParams{j}  = net.kernelSize(j) * net.kernelSize(j) * net.filter(j) * net.filter(j+1) + net.filter(j+1);
%                         idxb{j} = colon(1, net.filter(j+1))';   
                    else                     
                        numParams{j}  = net.kernelSize(j) *...
                            net.kernelSize(j) * net.filter(j) *...
                            net.filter(j+1) + net.filter(j+1);
                        idxb{j} = colon(1, net.filter(j+1))';                       
                    end
                    idxw{j} = colon(1, net.kernelSize(j) *...
                        net.kernelSize(j) * net.filter(j) *...
                        net.filter(j+1))';
                elseif layer(j+1) == net.layerEncoder.tconv % Transposed conv layer
                    if net.paddingType(j) == 1
                        imgWloop = double(net.stride(j)) *...
                            (double(net.imgW(j)) - 1) +...
                            double(net.kernelSize(j)) -...
                            2 * double(net.padding(j));
                        imgHloop = double(net.stride(j)) *...
                            (double(net.imgH(j)) - 1) +...
                            double(net.kernelSize(j)) -...
                            2 * double(net.padding(j));
                    elseif net.paddingType(j) == 2
                        imgWloop = double(net.stride(j)) *...
                            (double(net.imgW(j)) - 1) +...
                            double(net.kernelSize(j)) -...
                            double(net.padding(j));
                        imgHloop = double(net.stride(j)) *...
                            (double(net.imgH(j)) - 1) +...
                            double(net.kernelSize(j)) -...
                            double(net.padding(j));
                    end
                    if floor(imgWloop)~=imgWloop ...
                            || floor(imgHloop)~=imgHloop
                        error('The hyperparameters for conv. layer are invalid')
                    else
                        net.imgW(j+1) = imgWloop;
                        net.imgH(j+1) = imgHloop;
                    end
                    if j < numLayers-2 ...
                            && (layer(j+2) == net.layerEncoder.ln ...
                            || layer(j+2) == net.layerEncoder.bn) 
                        numParams{j} = net.kernelSize(j) * net.kernelSize(j)...
                            * net.filter(j) * net.filter(j+1);                    
                    else                     
                        numParams{j}  = net.kernelSize(j) *...
                            net.kernelSize(j) * net.filter(j) *...
                            net.filter(j+1) + net.filter(j+1);
                        idxb{j} = colon(1, net.filter(j+1))';                       
                    end
                    idxw{j} = colon(1, net.kernelSize(j) * net.kernelSize(j)...
                        * net.filter(j) * net.filter(j+1))';
                elseif layer(j+1) == net.layerEncoder.mp ...
                        || layer(j+1) == net.layerEncoder.ap % pooling layer
                    if net.paddingType(j) == 1
                        imgWloop = (double(net.imgW(j)) -...
                            double(net.kernelSize(j)) +...
                            2 * double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                        imgHloop = (double(net.imgH(j)) -...
                            double(net.kernelSize(j)) +...
                            2 * double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                    else
                        imgWloop = (double(net.imgW(j)) -...
                            double(net.kernelSize(j)) +...
                            double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                        imgHloop = (double(net.imgH(j)) -...
                            double(net.kernelSize(j)) +...
                            double(net.padding(j))) /...
                            double(net.stride(j)) + 1;
                    end
                    if floor(imgWloop)~=imgWloop ...
                            || floor(imgHloop)~=imgHloop
                        error('The hyperparameters for pooling. layer are invalid')
                    else
                        net.imgW(j+1) = imgWloop;
                        net.imgH(j+1) = imgHloop;
                    end
                elseif layer(j+1) == net.layerEncoder.ln ...
                        || layer(j+1) == net.layerEncoder.bn % Layer and batch normalization
                    net.imgW(j+1) = net.imgW(j);
                    net.imgH(j+1) = net.imgH(j);
                    if layer(j) == net.layerEncoder.fc
                        numParams{j} = 2 * nodes(j);
                    else
                        numParams{j} = 2 * net.filter(j);
                    end
                    idxb{j}  = colon(1, numParams{j} / 2)';
                    idxw{j}  = colon(1, numParams{j} / 2)';
                end
                if ~isempty(idxw{j})
                    numParamsPerlayer(1, j) = length(idxw{j});
                    paramUpdateIdx(1, j)    = 1;
                end
                if ~isempty(idxb{j})
                    numParamsPerlayer(2, j) = length(idxb{j});
                    paramUpdateIdx(2, j)    = 1;
                end  
                if net.xsc(j+1) ~= 0 ...
                        && (net.filter(net.xsc(j+1)) ~= net.filter(j+1) ...
                        ||net.imgW(net.xsc(j+1)) ~= net.imgW(j+1))
                    idxXsc = net.xsc(j+1);
                    idxwXsc{idxXsc} = colon(1, net.filter(idxXsc) *...
                        net.filter(j+1))';
                    idxbXsc{idxXsc} = colon(1, net.filter(j+1))';
                    numParams{j} = numParams{j} + net.filter(idxXsc) *...
                        net.filter(j+1) + net.filter(j+1);
                    if ~isempty(idxwXsc{idxXsc})
                        numParamsPerlayer(3, idxXsc) = length(idxwXsc{idxXsc});
                        paramUpdateIdx(3, idxXsc) = 1;
                    end
                    if ~isempty(idxbXsc{idxXsc})
                        numParamsPerlayer(4, idxXsc) = length(idxbXsc{idxXsc});
                        paramUpdateIdx(4, idxXsc) = 1;
                    end
                end 
                if ~isempty(numParams{j})
                    totalNumParams = totalNumParams +...
                        cast(numParams{j}, class(totalNumParams));
                end                            
            end   
            numParamsPerlayer_2 = [zeros(4, 1, 'like', numParamsPerlayer),...
                numParamsPerlayer];
            net.nodes   = nodes;
            net.idxb    = idxb;
            net.idxw    = idxw;
            net.idxwXsc = idxwXsc;
            net.idxbXsc = idxbXsc;
            net.totalNumParams = totalNumParams;
            net.numParamsPerlayer = numParamsPerlayer; % to be removed in the next version
            net.numParamsPerlayer_2 = cumsum(numParamsPerlayer_2, 2);
            net.similarIdx = []; % indices.getSimilarity(net); commented out since used by conv layers only
            net.paramUpdateIdx = paramUpdateIdx;
            
            % Transfer to gpu
%             if net.gpu==1
%                 net.numParamsPerlayer_2  = gpuArray(net.numParamsPerlayer_2);
%             end
        end
        function net = covariance(net)
            % Initialization
            batchSize   = net.batchSize;
            numLayers   = length(net.nodes);
            layer       = net.layer;
            % Indices for F*mwa
            idxFmwa     = cell(numLayers - 1, 2);
            idxFmwaXsc  = cell(numLayers - 1, 2);
            % Indices for F*Czwa
            idxFCzwa    = cell(numLayers - 1, 2);
            idxFCzwaXsc = cell(numLayers - 1, 2);
            % Indices for F*Cwz
            idxFCwz     = cell(numLayers - 1, 2);
            % Indices for updating hidden states between layers
            idxSzzUd    = cell(numLayers - 1, 1);
            idxSzzUdXsc = cell(numLayers - 1, 1);
            % Indices for the next hidden states used to updated w and b
            idxSwzUd    = cell(numLayers - 1, 1);
            % Padding
            paddingXsc  = zeros(1, numLayers);
            % Indices for the pooling layer
            idxPooling  = cell(numLayers - 1, 1);  
            
            idxFmwaXscT  = cell(numLayers - 1, 2);
            idxFCzwaXscT = cell(numLayers - 1, 2);
            idxSzzUdXscT = cell(numLayers - 1, 1);
            idxPoolingT  = cell(numLayers - 1, 1);
            idxFmwaT     = cell(numLayers - 1, 2);
            idxFCzwaT    = cell(numLayers - 1, 2);
            idxSzzUdT    = cell(numLayers - 1, 1);
            idxSwzUdT    = cell(numLayers - 1, 1);
            idxFCwzT     = cell(numLayers - 1, 2);
            for j = 1:numLayers - 1 
                % Conv. layer 
                if layer(j+1) == net.layerEncoder.conv ...
                        && net.similarIdx(j) == j                      
                    [idxFmwa{j, 2}, idxFCzwa{j, 1}, ~,...
                        idxSzzUd{j}] = indices.conv(net.kernelSize(j),...
                        net.stride(j), net.imgW(j), net.imgH(j), net.filter(j),...
                        net.imgW(j+1), net.imgH(j+1), net.filter(j+1),...
                        net.padding(j), net.paddingType(j), net.idxw{j},...
                        batchSize, net.dtype); 
                    idxFCzwa{j, 1} = idxFCzwa{j, 1}';                   
                    if net.cuda
                        n               = size(idxFCzwa{j, 1}, 1);  
                        idxFmwaT{j, 2}  = idxFmwa{j, 2}(1 : net.kernelSize(j)...
                            * net.kernelSize(j), 1 : net.imgW(j+1)...
                            * net.imgH(j+1));
                        idxFCzwaT{j, 1} = idxFCzwa{j, 1}(1 :(n /...
                            double(net.filter(j+1))), 1 : net.imgW(j) *...
                            net.imgH(j));
                        idxSzzUdT{j}    = idxSzzUd{j}(1 : net.imgW(j) *...
                            net.imgH(j), 1 :(n / double(net.filter(j+1))));
                    end
                    if j==1&&~net.convariateEstm
                        idxSzzUd{j} = 0;
                        idxFCzwa{j, 1} = 0;
                    end
                % Transposed cpnv. layer
                elseif layer(j+1) == net.layerEncoder.tconv 
                    [~, idxFCzwa_1_conv, idxFCzwa_2_conv,...
                        idxSzzUd_conv] = indices.conv(net.kernelSize(j),...
                        net.stride(j), net.imgW(j+1), net.imgH(j+1),...
                        net.filter(j+1), net.imgW(j), net.imgH(j),...
                        net.filter(j), net.padding(j), net.paddingType(j),...
                        net.idxw{j}, batchSize, net.dtype); 
                    idxFmwa{j, 1} = idxFCzwa_1_conv';
                    idxFmwa{j, 2} = idxSzzUd_conv';
                    [idxFCwz{j, 2}, idxSwzUd{j}, idxFCzwa{j, 1},...
                        idxSzzUd{j}] = indices.tconv(idxFCzwa_1_conv,...
                        idxFCzwa_2_conv, idxSzzUd_conv, net.kernelSize(j),...
                        net.imgW(j+1), net.imgH(j+1), net.filter(j+1),...
                        net.imgW(j), net.imgH(j), net.filter(j), batchSize);
                    if net.cuda
                        n               = size(idxFmwa{j, 1});
                        n2              = size(idxFCwz{j, 2});   
                        idxFmwaT{j, 1}  = idxFmwa{j, 1}(1 : n(1) /...
                            net.filter(j), 1 : net.imgW(j+1) * net.imgH(j+1));
                        idxFmwaT{j, 2}  = idxFmwa{j, 2}(1 : n(1) /...
                            net.filter(j), 1 : net.imgW(j+1) * net.imgH(j+1));
                        idxSzzUdT{j}    = idxSzzUd{j}(1 : net.kernelSize(j) *...
                            net.kernelSize(j), 1 : net.imgW(j) * net.imgH(j));
                        idxFCzwaT{j, 1} = idxFCzwa{j, 1}(1 : net.kernelSize(j) *...
                            net.kernelSize(j), 1 : net.imgW(j) * net.imgH(j));                                          
                        idxSwzUdT{j}    = idxSwzUd{j}(1 : n2(1) /...
                            batchSize, 1 : net.kernelSize(j) * net.kernelSize(j));
                        idxFCwzT{j, 2}  = idxFCwz{j, 2}(1 : n2(1) /...
                            batchSize, 1 : net.kernelSize(j) * net.kernelSize(j));
                    end
                elseif layer(j+1) == net.layerEncoder.mp ...
                        || layer(j+1) == net.layerEncoder.ap 
                    if net.kernelSize(j) == net.stride(j) ...
                            ||(net.kernelSize(j)==net.imgW(j) ...
                            && net.stride(j)==1)
                        overlap = 0;
                    else
                        overlap = 1;
                    end
                    % Pooling indices
                    paddingIdx        = net.imgW(j) * net.imgH(j) *...
                        net.filter(j) * batchSize + 1;
                    [img,...
                        paddingImg] = indices.imageConstruction(net.imgW(j),...
                        net.imgH(j), net.padding(j), paddingIdx,...
                        net.paddingType(j));
                    idxPooling_1      = indices.receptiveField(img,...
                        net.kernelSize(j), net.stride(j),  net.imgW(j+1),...
                        net.imgH(j+1)); 
                    % Copie for each filter
                    filterType        = 2;
                    [idxPooling_1,...
                        rawIdx] = indices.eachFilter(idxPooling_1,...
                        paddingImg, net.imgW(j), net.imgH(j), net.filter(j),...
                        net.padding(j), paddingIdx, filterType);                    
                    if layer(j+1) == net.layerEncoder.ap ...
                            && overlap == 1
                        baseIdx       = indices.idemSort(idxPooling_1,...
                            net.padding(j));
                        baseIdx       = repmat(baseIdx, [batchSize, 1]);
                        baseIdxM      = cast(baseIdx, class(net.imgW(j)));
                        baseIdxM      = baseIdxM .* (colon(1, net.imgW(j) *...
                            net.imgH(j) * net.filter(j) * batchSize)');
                        idxFCzwa_2    = baseIdxM(baseIdx);                        
                    end
                    % Copie for each batch
                    typeBatch         = 2;
                    idxPooling_1      = indices.eachBatch(idxPooling_1,...
                        rawIdx, net.imgW(j), net.imgH(j), net.filter(j),...
                        batchSize, typeBatch);                    
                    if layer(j+1) == net.layerEncoder.ap ...
                            && overlap == 1
                        [idxFCzzRef,...
                            idxNoPadding] = indices.refSort(idxFCzwa_2,...
                            idxPooling_1, net.padding(j), paddingIdx);
                        idxSzzUd{j} = reshape(repmat(colon(1, net.imgW(j+1) *...
                            net.imgH(j+1) * net.filter(j+1) * batchSize),...
                            [net.kernelSize(j) * net.kernelSize(j), 1]),...
                            [net.imgW(j+1) * net.imgH(j+1) * net.filter(j+1)...
                            .* net.kernelSize(j) * net.kernelSize(j)...
                            * batchSize, 1]);
                        if net.padding(j) ~= 0  
                            idxSzzUd{j}(idxNoPadding) = [];
                        end
                        idxSzzUd{j}        = idxSzzUd{j}(idxFCzzRef);
                        baseIdxM(baseIdx)  = idxSzzUd{j};
                        baseIdxM(~baseIdx) = net.imgW(j+1) * net.imgH(j+1)...
                            * net.filter(j+1) * batchSize + 1;
                        idxSzzUd{j}        = cast(baseIdxM, class(idxSzzUd{j}));
                        if net.cuda
                            idxSzzUdT{j} = idxSzzUd{j}(1 : net.imgW(j) * net.imgH(j), :);
                        end
                    else %Attention !!! need to fix or remove this
                        if net.cuda
                            idxSzzUdT{j} = gpuArray(single(1));
                        end
                    end                  
                    idxPooling{j} = idxPooling_1;
                    if net.cuda
                        idxPoolingT{j} = idxPooling{j}(1 : net.imgW(j+1)...
                            * net.imgH(j+1), 1 : net.kernelSize(j) * net.kernelSize(j));
                    end
                    clearvars idxNoPadding rawIdx idxPooling_1
                end 
                % Indices for shortcut in residual network
                if net.xsc(j) ~= 0 ...
                        && (net.filter(net.xsc(j)) ~= net.filter(j) ...
                        || net.imgW(net.xsc(j)) ~= net.imgW(j))
                    idxXsc             = net.xsc(j);
                    paddingXsc(idxXsc) = 0; 
                    kernelSizeXsc      = 1;
                    paddingTypeXsc     = 1;
                    [idxFmwaXsc{idxXsc, 2}, idxFCzwaXsc{idxXsc, 1},...
                        idxFCzwaXsc{idxXsc, 2}, idxSzzUdXsc{idxXsc, 1}...
                        ] = indices.conv(kernelSizeXsc, net.stride(idxXsc),...
                        net.imgW(idxXsc), net.imgH(idxXsc), net.filter(idxXsc),...
                        net.imgW(j), net.imgH(j), net.filter(j),...
                        paddingXsc(idxXsc), paddingTypeXsc,...
                        net.idxwXsc{idxXsc}, batchSize, net.dtype);   
%                     idxSzzUdXsc{idxXsc, 1} = idxSzzUdXsc{idxXsc, 1}(:, 1);
                    idxFCzwaXsc{idxXsc, 2} = idxFCzwaXsc{idxXsc, 2}(:,1);
                    if net.gpu
                        idxFmwaXsc{idxXsc, 1}  = gpuArray(single(idxFmwaXsc{idxXsc, 1}));
                        idxFmwaXsc{idxXsc, 2}  = gpuArray(single(idxFmwaXsc{idxXsc, 2}));
                        
                        idxFCzwaXsc{idxXsc, 1} = gpuArray(single(idxFCzwaXsc{idxXsc, 1}));
                        idxFCzwaXsc{idxXsc, 2} = gpuArray(single(idxFCzwaXsc{idxXsc, 2}));
                        idxSzzUdXsc{idxXsc}    = gpuArray(single(idxSzzUdXsc{idxXsc}));
                    end
                    if net.cuda %!!! attention
                        idxFmwaXscT{idxXsc, 2}  = idxFmwaXsc{idxXsc, 2}(1,...
                            1 : net.imgW(j+1) * net.imgH(j+1));
                        idxFCzwaXscT{idxXsc, 1} = idxFCzwaXsc{idxXsc, 1}';
                        n                       = size(idxFCzwaXscT{idxXsc, 1}, 1);
                        q2                      = size(idxFCzwaXscT{idxXsc, 1}, 2);
                        wh2                     = sqrt(q2 / double(net.filter(idxXsc)));
                        idxFCzwaXscT{idxXsc, 1} = idxFCzwaXscT{idxXsc, 1}(1 : (n / double(net.filter(j+1))), 1 : wh2 * wh2);
                        idxSzzUdXscT{idxXsc, 1} = idxSzzUdXsc{idxXsc, 1}(1 : wh2 * wh2, 1 : (n / double(net.filter(j+1))));
                        idxFCzwaXscT{idxXsc, 2} = idxFCzwaXsc{idxXsc, 2}(1 : wh2 * wh2, :);
                    end
                end
                if net.gpu
                    idxFmwa{j, 1}  = gpuArray(single(idxFmwa{j, 1}));
                    idxFmwa{j, 2}  = gpuArray(single(idxFmwa{j, 2}));                                     
                    idxFCzwa{j, 1} = gpuArray(single(idxFCzwa{j, 1}));
                    idxFCzwa{j, 2} = gpuArray(single(idxFCzwa{j, 2}));
                    idxSzzUd{j}    = gpuArray(single(idxSzzUd{j}));
                    idxPooling{j}  = gpuArray(single(idxPooling{j}));  
                    idxFCwz{j}     = gpuArray(idxFCwz{j});
                    idxSwzUd{j}    = gpuArray(idxSwzUd{j});
                    
                    idxFmwaT{j, 1}  = gpuArray(single(idxFmwaT{j, 1}));
                    idxFmwaT{j, 2}  = gpuArray(single(idxFmwaT{j, 2}));                                     
                    idxFCzwaT{j, 1} = gpuArray(single(idxFCzwaT{j, 1}));
%                     idxFCzwaT{j, 2} = gpuArray(single(idxFCzwaT{j, 2}));
                    idxSzzUdT{j}    = gpuArray(single(idxSzzUdT{j}));
                    idxPoolingT{j}  = gpuArray(single(idxPoolingT{j})); 
                    idxFCwzT{j, 2}  = gpuArray(single(idxFCwzT{j, 2}));
                    idxSwzUdT{j}    = gpuArray(single(idxSwzUdT{j}));
                end
            end
            % Outputs
            net.idxFmwa    = idxFmwa;
            net.idxFmwaXsc = idxFmwaXsc;
            net.idxFCzwa   = idxFCzwa;
            net.idxFCzwaXsc= idxFCzwaXsc;
            net.idxFCwz    = idxFCwz;
            net.idxSzzUd   = idxSzzUd;          
            net.idxSzzUdXsc= idxSzzUdXsc;
            net.idxSwzUd   = idxSwzUd;
            net.idxPooling = idxPooling;
            net.paddingXsc = paddingXsc; 
            
            net.idxFmwaT     = idxFmwaT;
            net.idxFCzwaT    = idxFCzwaT;
            net.idxSzzUdT    = idxSzzUdT;
            net.idxPoolingT  = idxPoolingT;
            
            net.idxFmwaXscT  = idxFmwaXscT;
            net.idxFCzwaXscT = idxFCzwaXscT;
            net.idxSzzUdXscT = idxSzzUdXscT;
            net.idxSwzUdT    = idxSwzUdT;
            net.idxFCwzT     = idxFCwzT;
        end
    end
end