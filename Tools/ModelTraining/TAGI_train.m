function [best_theta, best_LL, terminate] = TAGI_train(net, xtrain, ytrain, S_ytrain, xvalid, yvalid, seed, load_theta)
    if net.gpu
        xtrain  = gpuArray(single(xtrain));
        ytrain  = gpuArray(single(ytrain));
        S_ytrain  = gpuArray(single(S_ytrain));
        xvalid   = gpuArray(single(xvalid));
        yvalid   = gpuArray(single(yvalid));    
    end
    
    % Set the seed
    if seed
        rng(seed);
        net.seed = seed;
    else
        rng('shuffle');
    end

    maxEpoch   = net.maxEpoch;

    % Values used to find the optimal number of epochs
    best_LL = -inf;
    patience = 0;
    patience_threshold = 1;

    % Train net
    net.trainMode = 1;
    [net, states, maxIdx, netInfo] = network.initialization(net);
    normStat = tagi.createInitNormStat(net);

    % Validation net
    netT              = net;
    netT.trainMode    = 0;
    [netT, statesT, maxIdxT] = network.initialization(netT); 
    normStatT = tagi.createInitNormStat(netT);
    
    % Initalize weights and biases
    if ~isempty(load_theta)
        theta = load_theta;
    else
        theta = tagi.initializeWeightBias(net);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Training
    stop  = 0;
    epoch = 0;
    reinit = 0;
    best_epoch = 0;
    terminate = false;
    fprintf('Training BNN using TAGI...\n');
    while best_LL == -inf
        if stop
            if reinit == 0
                fprintf('\nBNN did not converge... re-initializing weights and biases...');
            else
                fprintf('...');
            end
            theta = tagi.initializeWeightBias(net);
            epoch = 0;
            stop = 0;
            patience = 0;
            reinit = reinit + 1;
        end
        
        if reinit > 2
            fprintf('\nRe-initializing of weights and biases exceeded 3 attempts... Terminating...\n');
            if ~isempty(load_theta)
                fprintf('BNN failed to learn... Returning the previously loaded weights and biases...\n');
                best_theta = load_theta;
                terminate = true;
            else
                fprintf('BNN failed to learn... Try different network configuration, data scaling, or weight and bias initialization ...\n');
                best_theta = [];
            end
            break
        end
        
        while ~stop
            if epoch >1
                idxtrain = randperm(size(ytrain, 1));
                ytrain   = ytrain(idxtrain, :);
                xtrain   = xtrain(idxtrain, :);
            end

            [theta, normStat, mytrain, Sytrain,...
                ~] = network.regression(net, theta, normStat, states,...
                maxIdx, xtrain, ytrain, S_ytrain);

            if net.learnSv
                [mv2a, ~, ~] = act.expFun(mytrain(:,2), Sytrain(:,2), net.gpu);
                Sytrain(:,1) = Sytrain(:,1) +  mv2a;
            end
            
            if any(isnan(mytrain), 'all')
                fprintf('\nWarn! NaN values in mytrain\n');
            end
            if any(isnan(Sytrain), 'all')
                fprintf('\nWarn! NaN values in Sytrain\n');
            end
            
            % Validating
            [~, ~, myvalid, Syvalid] = network.regression(netT, theta, ...
            normStatT, statesT, maxIdxT, xvalid, []);
        
            if netT.learnSv
                [mv2a, ~, ~] = act.expFun(myvalid(:,2), Syvalid(:,2), net.gpu);
                Syvalid(:,1) = Syvalid(:,1) +  mv2a;
            end
        
            LL_val = mt.loglik(yvalid, myvalid(:,1), Syvalid(:,1));
            if isnan(LL_val)
                LL_val = -inf;
            end

            if LL_val > best_LL
                best_LL = LL_val;
                best_theta = theta;
                best_epoch = epoch;
                patience = 0;
                reinit = 0;
            else
                patience = patience + 1;
            end
            
            epoch = epoch + 1;
            
            if epoch > maxEpoch
                fprintf('Epoch > %d, max epoch exceeded. Training terminated... best epoch: %d \n', maxEpoch, best_epoch)
                stop = true;
            elseif patience == patience_threshold
                fprintf('Early stopping... Training terminated at epoch %d, best epoch: %d \n', epoch-1, best_epoch)
                stop = true;
            end
        end
    end
    fprintf('************************************\n')
end