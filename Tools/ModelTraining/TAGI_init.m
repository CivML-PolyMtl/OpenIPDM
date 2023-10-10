function net = TAGI_init(modelName, maxEpoch, batchSize, x_dim, y_dim, mw_gain, Sw_gain, mb_gain, Sb_gain, h_nodes, h_layer_act, o_layer_act)
    % Set NN attributes
    net.task           = 'regression';
    net.modelName      = modelName;
    net.saveModel      = 0;
    net.maxEpoch       = maxEpoch;
    % GPU
    net.gpu            = false;
    net.cuda           = false;
    net.numDevices     = 1;
    % Data type object half or double precision
    net.dtype          = 'double';
    % Number of input covariates
    net.nx             = x_dim; 
    % Number of output responses
    net.nl             = y_dim;
    net.nv2            = y_dim;
    net.ny             = 1*y_dim; 
    % Batch size 
    net.batchSize      = batchSize; 
    net.repBatchSize   = 1;

    % Noise parameter
    net.sv = 0.1;
    net.learnSv = 1;
    net.split_v2_act = 0;
    if net.learnSv
        net.ny = 2*y_dim;
        net.sv = 0;
        net.noiseType = 'hete';
        net.gain_v2 = false;
%         net.noiseType = 'full';
    end

    if strcmp(net.noiseType, 'full')
        net.nLchol         = (net.nv2^2 + net.nv2)/2;
        net.ny             = net.nl + net.nLchol; 
        net.gain_v2 = false;
    else
        net.nLchol = y_dim;
        net.gain_v2 = false;
    end

    % Build the network architecture
    % Layer| 1: FC; 2:conv; 3: max pooling; 4: avg pooling; 5: layerNorm; 6: batchNorm 
    % Activation: 1:tanh; 2:sigm; 3:cdf; 4:relu; 5:softplus
    net.layer          = [1         1  1];
    net.nodes          = [net.nx    h_nodes net.ny]; 
    net.actFunIdx      = [0         h_layer_act o_layer_act];

    % Parameter initialization method
    net.initParamType  = 'He';
    numLayers = length(net.layer);
    if mw_gain
        net.gainMw = mw_gain.*ones(1, numLayers-1);
    end
    if Sw_gain
        net.gainSw = Sw_gain.*ones(1, numLayers-1);
    end
    if mb_gain
        net.gainMb = mb_gain*ones(1, numLayers-1);
    end
    if Sb_gain
        net.gainSb = Sb_gain*ones(1, numLayers-1);
    end
%     net.gainMv2=1e-2;
%     net.gainSv2=1e-3;
end