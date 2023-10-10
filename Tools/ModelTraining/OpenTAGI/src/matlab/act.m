%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File:         act
% Description:  Activation function
% Authors:      Luong-Ha Nguyen & James-A. Goulet & Ali Fakhri 
% Created:      November 12, 2019
% Updated:      September 15, 2023 | Removed functions not used by OpenIPDM
% Contact:      luongha.nguyen@gmail.com & james.goulet@polymtl.ca & afakhri@pm.me
% Copyright (c) 2021 Luong-Ha Nguyen & James-A. Goulet. All rights reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 classdef act   
    methods (Static)
        function [m, S, J] = meanVar(z, mz, Sz, funIdx, bound, B, rB, gpu)         
            if funIdx == 1 % tanh
                if gpu
                    dtanhf = @(x) (1 - tanh(x).^2);
                    m = tanh(mz) * bound;
                    J = arrayfun(dtanhf, mz);
                    J = J .* bound;
                else
                    dtanhf = @(x) 1 - tanh(x).^2;
                    m = dtanhf(mz) .* (z - mz) + tanh(mz);
                    J = dtanhf(z);  
                end
            elseif funIdx == 2 % sigmoid
                if gpu
                    sigmoid_mz  = exp(-mz);
                    sigmoid_mz  = bsxfun(@plus, 1, sigmoid_mz);
                    sigmoid_mz  = bsxfun(@rdivide, 1, sigmoid_mz);                  
                    dsigmoid_mz = bsxfun(@minus, 1, sigmoid_mz);
                    dsigmoid_mz = bsxfun(@times, sigmoid_mz, dsigmoid_mz);                                     
                    m  = sigmoid_mz;                     
                    J  = dsigmoid_mz;  
                else
                    sigmoid  = @(x) 1 ./ (1 + exp(-x));
                    dsigmoid = @(x) sigmoid(x) .* (1 - sigmoid(x));
                    m = sigmoid(mz);
                    J = dsigmoid(z);
                end
            elseif funIdx == 3 % cdf
                if gpu
                    m = normcdf(mz);
                    J = normpdf(mz);
                else
                    m  = normpdf(mz) .* (z - mz) + normcdf(mz);
                    J  = normpdf(z);
                end
            elseif funIdx == 4 % relu
                if gpu
                    J = mz > 0;
                    J   = cast(J, 'like', mz);
                    m   = bsxfun(@times, z, J);
                else
                    m   = max(0, mz);
                    J   = single(z > 0);
                end            
            elseif funIdx == 5 % softplus
                if gpu
                    alpha = 2;
                    k = alpha * mz < 1000;
                    e = bsxfun(@plus, 1, exp(alpha * mz .* k));
                    m = (log(e) + mz .* (1 - k)) / alpha;
                    J = k .* bsxfun(@rdivide, exp(alpha * mz .* k), e) ...
                        + (1 - k);
                else
                    m = log(1 + exp(mz));
                    J = 1 ./ (1 + exp(-mz));
                end
            elseif funIdx == 6 % leaky relu
                alpha = cast(0.2, 'like', mz);
                if gpu                   
                    idx = mz > 0;
                    J   = cast(idx, 'like', mz);                   
                    m   = bsxfun(@times, z, J);  
                    J(~idx) = alpha;
                    m(~idx) = alpha * z(~idx);
                else
                    idx = mz > 0;
                    m   = max(0, mz);
                    J   = single(z > 0);
                    J(~idx) = alpha;
                    m(~idx) = alpha * z(~idx);
                end 
            elseif funIdx == 7 % exponential relu
                 alpha = cast(0.001, 'like', mz);
                if gpu
                    idx = mz > 0;
                    m   = mz;
                    m(~idx ) = alpha * (exp(mz(~idx) + 0.5*Sz(~idx)) - 1);
                    J = cast(idx, 'like', mz);  
                    J(~idx) = alpha*exp(mz(~idx) + 0.5 * Sz(~idx));                          
                else
                    idx = mz > 0;
                    m   = mz;
                    m(~idx ) = alpha * (exp(mz(~idx) + 0.5 * Sz(~idx)) - 1);
                    J = cast(idx, 'like', mz);  
                    J(~idx) = alpha * exp(mz(~idx) + 0.5 * Sz(~idx));
                end 
            elseif funIdx == 8
                if gpu
                    m = 1*sin(mz);
                    J = 1*cos(mz);
                else
                    m = sin(mz);
                    J = cos(mz);
                end
            elseif funIdx == 9
                alpha = 2;
                if gpu
                    dtanhf = @(x) 1 - tanh(x).^2;
                    m = alpha .* (tanh(mz) + 1);
                    J = arrayfun(dtanhf, mz);
                    J = alpha .* J;
                else
                    dtanhf = @(x) 1-tanh(x).^2;
                    m = alpha .* (tanh(mz) + 1);
                    J = arrayfun(dtanhf, mz);
                    J = alpha .* J;
                end 
            elseif funIdx == 10 % softmax
                ny = length(mz)/(B * rB);
                mz = reshape(mz, [ny, B * rB]);
                if gpu
                    maxMz   = max(mz);
                    mzShift = bsxfun(@minus, mz, maxMz);
                    expMz   = exp(mzShift);
                    m       = bsxfun(@rdivide, expMz, sum(expMz));
                    m       = m(:);
                    fun     = @(x) (1 - x) .* x;
                    J       = arrayfun(fun, m);
                else
                    maxMz   = max(mz);
                    mzShift = bsxfun(@minus, mz, maxMz);
                    expMz   = exp(mzShift);
                    m       = bsxfun(@rdivide, expMz, sum(expMz));
                    fun     = @(x) (1 - x) .* x;
                    J       = arrayfun(fun, m);
                end
            elseif funIdx == 11 % relu
                alpha = cast(-0.2, 'like', mz);
                if gpu                   
                    idx = mz < 0;
                    J   = cast(idx, 'like', mz);                   
                    m   = bsxfun(@times, z, J);  
                    J(~idx) = alpha;
                    m(~idx) = alpha * z(~idx);
                else
                    idx = mz < 0;
                    m   = min(0, mz);
                    J   = -single(z < 0);
                    J(~idx) = alpha;
                    m(~idx) = alpha * z(~idx);
                end 
            else
                m = mz;
                J = ones(size(mz), 'like', mz);
            end 
            if gpu
                fun = @(x, y) (x.^2) .* y;            
                S   = arrayfun(fun, J, Sz);
            else
                S = J .* Sz .* J;
            end
            if funIdx == 7
                S(~idx) = alpha^2 * exp(2 * mz(~idx) + Sz(~idx)) ...
                    .* (exp(Sz(~idx)) - 1);
            end
        end
        function [ma, Sa, Cza] = expFun(mz, Sz, gpu)
            if gpu == 1
                ma  = exp(mz + 0.5 * Sz);
                Sa  = exp(2 * mz + Sz) .* (exp(Sz) - 1);
                Cza = Sz .* exp(mz + 0.5 * Sz);
            else
                ma  = exp(mz + 0.5 * Sz);
                Sa  = exp(2 * mz + Sz) .* (exp(Sz) - 1);
                Cza = Sz .* exp(mz + 0.5 * Sz);
            end
        end
    end
end