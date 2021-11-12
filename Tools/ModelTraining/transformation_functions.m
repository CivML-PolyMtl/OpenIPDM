%% Transformation functions
% Created by: James-A. Goulet
% November 22, 2016
%
%% INPUTS:
% bounds:      NxD cell containing lower and upper bound for each parameter
%              {[b11 b12],[b21,b22],...[bD1,bD2]}
% 
%% OUTPUTS:
% TR_fct:      original to transformed space function
% TR_fct_inv:  transformed to original space function
% dTR_fct:     1st derivative of TR_fct
% dTR_fct_inv: 1st derivative of TR_fct_inv
%%
function [TR_fct,TR_fct_inv,dTR_fct,dTR_fct_inv]=transformation_functions(bounds)

%% Analysis parameters
D=size(bounds,2); %number of parameters

%% String initialization
TR_fct_string=[]; %function string evaluation
TR_fct_inv_string=[]; %function string evaluation
dTR_fct_string=[]; %function string evaluation
dTR_fct_inv_string=[]; %function string evaluation

%% Function builder
for d=1:D
    d_s=num2str(d);
    if isempty(bounds{d})|| all(isinf(bounds{d})) %Case 1: x^{TR}=x (no trandformation)
        TR_fct_string=[TR_fct_string 'x(:,' d_s ') ']; 
        TR_fct_inv_string=[TR_fct_inv_string 'x(:,' d_s ') '];
        dTR_fct_string=[dTR_fct_string '1 ']; 
        dTR_fct_inv_string=[dTR_fct_inv_string '1 ']; 
    elseif bounds{d}(1)==0&&isinf(bounds{d}(2))   %Case 2: x^{TR}=log(x) (logarithmic transformation)
        TR_fct_string=[TR_fct_string 'log(x(:,' d_s ')) ']; 
        TR_fct_inv_string=[TR_fct_inv_string 'exp(x(:,' d_s ')) '];
        dTR_fct_string=[dTR_fct_string '1/x(:,' d_s ') ']; 
        dTR_fct_inv_string=[dTR_fct_inv_string 'exp(x(:,' d_s ')) ']; 
    elseif all(isfinite(bounds{d}))               %Case 3: x^{TR}=sigm(x) (sigmoid transdormation)
        a=num2str(bounds{d}(1));
        b=num2str(bounds{d}(2));
        TR_fct_string=[TR_fct_string '-log(((' b '-' a ')./(x(:,' d_s ')-' a '))-1) ']; 
        TR_fct_inv_string=[TR_fct_inv_string '((' b '-' a ')./(1+exp(-x(:,' d_s '))))+' a ' '];
        %dTR_fct_string=[dTR_fct_string '((' b '-' a ')./(1+exp(x(:,' d_s '))).^2)*exp(x(:,' d_s ')) ']; 
        %dTR_fct_inv_string=[dTR_fct_inv_string '(' a '-' b ')./((x(:,' d_s ')-' a ').*(x(:,' d_s ')-' b ')) '];
        dTR_fct_string=[dTR_fct_string '(' a '-' b ')./((x(:,' d_s ')-' a ').*(x(:,' d_s ')-' b ')) '];
        dTR_fct_inv_string=[dTR_fct_inv_string '((' b '-' a ')./(1+exp(-x(:,' d_s '))).^2)*exp(-x(:,' d_s ')) '];
    else
        error('Error: transformation_functions.m -> the bounds provided are not supported')
    end 
end

%% Transform strings into handle functions
eval(['TR_fct=@(x) [' TR_fct_string '];']);
eval(['TR_fct_inv=@(x) [' TR_fct_inv_string '];']);
eval(['dTR_fct=@(x) [' dTR_fct_string '];']);
eval(['dTR_fct_inv=@(x) [' dTR_fct_inv_string '];']);