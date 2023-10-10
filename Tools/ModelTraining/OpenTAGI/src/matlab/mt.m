%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         mt
% Description:  Metric for performance evaluation
% Authors:      Luong-Ha Nguyen & James-A. Goulet
% Created:      November 8, 2019
% Updated:      February, 2020
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef mt
    methods (Static)
        function e = computeError(y, ypred)
            e = mean(sqrt(mean((y-ypred).^2)));
        end
        function e  = errorRate(ytrue, ypred)
            idx_true     = ytrue;
            [~,idx_pred] = max(ypred);
            idx_pred     = idx_pred-1;
            e            = idx_true~=idx_pred;  
        end
        function LL = loglik(y, ypred, Vpred)
            d = size(y, 2);
            if d == 1
                LL = mean(-0.5*log(2*pi*Vpred) - (0.5*(y-ypred).^2)./Vpred);
            else
                LL = mean(-d/2*log(2*pi) - 0.5*log(prod(Vpred, 2)) - sum((0.5*(y-ypred).^2)./Vpred, 2)); 
            end
        end
        function weight     = probFromloglik(loglik)
            maxlogpdf       = max(loglik);
            w_1             = bsxfun(@minus,loglik,maxlogpdf);
            w_2             = log(sum(exp(w_1)));
            w_3             = bsxfun(@minus,w_1,w_2);
            weight          = exp(w_3);
        end
        function NLL = negLoglik(probPreds, labels)
            probPreds = gather(probPreds);
            labels    = gather(labels) + 1;
            numObs    = length(labels);
            classes   = colon(min(labels), max(labels));
            y1hot     = repmat(classes, [numObs, 1]);
            y1hot     = y1hot == labels;
            NLL       = -mean(log(probPreds(y1hot) + 1E-12));
        end
        function [ece, binLowers, avgConfsInBins] = expCaliError(probPreds, labels, nBins)
            probPreds = gather(probPreds);
            labels    = gather(labels) + 1;
            binBoundaries = linspace(0, 1, nBins + 1);
            binLowers = binBoundaries(1:nBins);
            binUppers = binBoundaries(2:nBins+1);
            
            % Initializaiton
            preds   = probPreds;
            targets = labels;
            [confidences, predictions] = max(preds, [], 2);
            accuracies = predictions == targets;
            eceloop = 0;
            for k = 1:nBins
                inBin = confidences > binLowers(k) & confidences <= binUppers(k);
                propInBin = mean(inBin);
                if propInBin > 0
                    accuracyInBin = mean(accuracies(inBin));
                    avgConfidenceInBin = mean(confidences(inBin));
                    delta = avgConfidenceInBin - accuracyInBin;
                    avgConfsInBins =  delta;
                    eceloop = eceloop + abs(delta) * propInBin;
                end
            end
            ece = eceloop;
        end
        function r = rocAucScore(probPreds, labels, numClasses)
            probPreds = gather(probPreds);
            labels    = gather(labels) + 1;
            rClasses  = zeros(numClasses, 1);
            for c = 1:numClasses
                pc = probPreds(:, c);
                [~,~,~,rClasses(c)] = perfcurve(labels, pc, c);
            end
            r = mean(rClasses);
        end
    end
end