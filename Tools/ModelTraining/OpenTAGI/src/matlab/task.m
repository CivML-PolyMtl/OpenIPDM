 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         task
% Description:  Run different tasks such as classification, regression, etc
% for the neural networks defined in the config file.
% Authors:      Luong-Ha Nguyen & James-A. Goulet & Ali Fakhri 
% Created:      July 02, 2020
% Updated:      September 15, 2023 | Removed functions not used by OpenIPDM
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca & afakhri@pm.me
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef task
    methods (Static)                                                                          
        % Regression
        function runRegressionFullCov(net, xtrain, Sxtrain, ytrain, xtest,...
                Sxtest, ytest)
            % Initialization          
            cdresults  = net.cd;
            modelName  = net.modelName;
            dataName   = net.dataName;
            saveModel  = net.saveModel;
            maxEpoch   = net.maxEpoch;
            svinit     = net.sv;
            
            % Train net
            net.trainMode = true;
            [net, states, maxIdx, netInfo] = network.initialization(net);
            normStat = tagi.createInitNormStat(net);
            
            % Test net
            netT              = net;
            netT.trainMode    = false;
            netT.batchSize    = 1;
            netT.repBatchSize = 1;
            [netT, statesT, maxIdxT] = network.initialization(netT); 
            normStatT = tagi.createInitNormStat(netT); 
            
            % Initalize weights and bias
            theta    = tagi.initializeWeightBias(net);            
            net.sv   = svinit;
  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Training
            stop  = 0;
            epoch = 0;
            tic
            while ~stop
                if epoch >1
                    idxtrain = randperm(size(ytrain, 1));
                    ytrain   = ytrain(idxtrain, :);
                    xtrain   = xtrain(idxtrain, :);
                    Sxtrain  = Sxtrain(idxtrain, :);
                    net.sv = 0.85*net.sv';
                    if net.sv<0.01
                        net.sv=0.01;
                    end
                end
                epoch = epoch + 1;
                [theta, normStat, mytrain, Sytrain,...
                    sv] = network.regressionFullCov(net, theta, normStat,...
                    states, maxIdx, xtrain, Sxtrain, ytrain);
                net.sv = sv;
                if epoch >= maxEpoch; break;end
            end
            toc
            Sytrain = Sytrain + net.sv.^2;
            figure('Position', [0 0 450 250]);
            [xtrain, idxS] = sort(xtrain);
            pl.regression(xtrain, mytrain(idxS), Sytrain(idxS), 'red', 'red', 1)
            hold on
            plot(xtrain, ytrain(idxS), 'k')
            hold off
            title('Training set')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Testing
            [~, ~, mytest, Sytest] = network.regressionFullCov(netT, theta,...
                normStatT, statesT, maxIdxT, xtest, Sxtest, []);
            Sytest = Sytest + net.sv.^2;
            figure('Position', [0 0 450 250]);
            [xtest, idxS] = sort(xtest);
            pl.regression(xtest, mytest(idxS), Sytest(idxS), 'red', 'red', 1)
            hold on
            plot(xtest, ytest(idxS), 'k')
            hold off
            title('Test set')
        end
        function runRegressionDiag(net, xtrain, Sxtrain, ytrain, xtest,...
                Sxtest, ytest)
            % Initialization          
            cdresults  = net.cd;
            modelName  = net.modelName;
            dataName   = net.dataName;
            saveModel  = net.saveModel;
            maxEpoch   = net.maxEpoch;
            svinit     = net.sv;
            
            % Train net
            net.trainMode = 1;
            [net, states, maxIdx, netInfo] = network.initialization(net);
            normStat = tagi.createInitNormStat(net);
            
            % Test net
            netT              = net;
            netT.trainMode    = 0;
            netT.batchSize    = 1;
            netT.repBatchSize = 1;
            [netT, statesT, maxIdxT] = network.initialization(netT); 
            normStatT = tagi.createInitNormStat(netT); 
            
            % Initalize weights and bias
            theta    = tagi.initializeWeightBias(net);            
            net.sv   = svinit;
  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Training
            stop  = 0;
            epoch = 0;
            tic
            while ~stop
                if epoch >1
                    idxtrain = randperm(size(ytrain, 1));
                    ytrain   = ytrain(idxtrain, :);
                    xtrain   = xtrain(idxtrain, :);
                    Sxtrain  = Sxtrain(idxtrain, :);
                    net.sv = 0.95*net.sv';
                    if net.sv<0.01
                        net.sv=0.01;
                    end
                end
                epoch = epoch + 1;
                [theta, normStat, mytrain, Sytrain,...
                    sv] = network.regression(net, theta, normStat, states,...
                    maxIdx, xtrain, ytrain);
                net.sv = sv;
                if epoch >= maxEpoch; break;end
            end
            toc
            Sytrain = Sytrain + net.sv.^2;
            figure('Position', [0 0 450 250]);
            [xtrain, idxS] = sort(xtrain);
            pl.regression(xtrain, mytrain(idxS), Sytrain(idxS), 'red', 'red', 1)
            hold on
            plot(xtrain, ytrain(idxS), 'k')
            hold off
            title('Training set')
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Testing
            [~, ~, mytest, Sytest] = network.regression(netT, theta, ...
                normStatT, statesT, maxIdxT, xtest, []);
            Sytest = Sytest + net.sv.^2;
            figure('Position', [0 0 450 250]);
            [xtest, idxS] = sort(xtest);
            pl.regression(xtest, mytest(idxS), Sytest(idxS), 'red', 'red', 1)
            hold on
            plot(xtest, ytest(idxS), 'k')
            hold off
            title('Test set')
        end      
                      
        % Save functions
        function saveRegressionNet(cd, modelName, dataName, theta, normStat,...
                metric, trainTime, netInfo, epoch)
            filename = [modelName, '_', 'E', num2str(epoch), '_', dataName];
            folder   = char([cd ,'/results/']);
            save([folder filename], 'theta', 'normStat', 'metric',...
                'trainTime', 'netInfo')
        end       
    end
end