%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         dp
% Description:  data processing
% Authors:      Luong-Ha Nguyen & James-A. Goulet & Ali Fakhri 
% Created:      November 8, 2019
% Updated:      September 15, 2023 | Removed functions not used by OpenIPDM
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca & afakhri@pm.me
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
classdef dp
    methods (Static)              
        % Normalization
        function [xntrain, yntrain, xntest, yntest, mxtrain, sxtrain,...
                mytrain, sytrain] = normalize(xtrain, ytrain, xtest, ytest)
            mxtrain = nanmean(xtrain);
            sxtrain = nanstd(xtrain);
            idx     = sxtrain==0;
            sxtrain(idx) = 1;
            mytrain = nanmean(ytrain);
            sytrain = nanstd(ytrain);
            xntrain = (xtrain - mxtrain)./sxtrain;
            yntrain = (ytrain - mytrain)./sytrain;
            xntest  = (xtest - mxtrain)./sxtrain;
            yntest  = ytest;
        end
        function [y, sy] = denormalize(yn, syn, myntrain, syntrain)
            y   = yn.*syntrain + myntrain;
            if ~isempty(syn)
                sy  = (syntrain.^2).*syn;
            else
                sy  = [];
            end
        end
        
        % Shared functions         
        function [xtrain, ytrain, xtest, ytest] = split(x, y, ratio)
            numObs      = size(x, 1);
            idxobs      = 1:numObs;
%             idxobs      = randperm(numObs);
            idxTrainEnd = round(ratio*numObs);
            idxTrain    = idxobs(1:idxTrainEnd);
            idxTest     = idxobs((idxTrainEnd+1):numObs);
            xtrain      = x(idxTrain, :);
            ytrain      = y(idxTrain, :);
            xtest       = x(idxTest, :);
            ytest       = y(idxTest, :);
        end
        function [trainIdx, testIdx] = indexSplit(numObs, ratio, dtype)
           idx = randperm(numObs);
           trainIdxEnd = round(numObs*ratio);
           trainIdx = idx(1:trainIdxEnd)';
           testIdx  = idx(trainIdxEnd+1:end)';
           if strcmp(dtype, 'single')
               trainIdx = int32(trainIdx);
               testIdx = int32(testIdx);
           end
        end
    end
end