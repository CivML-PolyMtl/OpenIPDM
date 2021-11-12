%% Newton-Raphson maximizer
% Created by: James-A. Goulet
% November 22, 2016
%
%% INPUTS:
% Mandatory inputs
% x0:    [Dx1] vector containing mean initial values for x
% fx:    @(x)  function for be maximized
%
% Optional inputs -> 'option name'
% output:  'original/transform'         %output results in original or
%                                        transformed space
%                                        Note: if the output is the
%                                        original space, the MLE is in the
%                                        original space too. if the output is the
%                                        transformed space, the MLE is in the
%                                        transformed space too.
% laplace: 'yes/'no'                    %perform laplace approximation
%                                        yes/no
% log_tranform:                         %perform NR in the log-transformed
%                                        space yes/no
% display: 'yes'/'no'                   %display results online
% nb_failed_iteration_limit             %maximal allowable number of failed iterations
% convergence_tol                       %convergence tolerance
% delta                                 %Relative delta x for numeric derivatives
% bounds: {[a1 bq],...,[aD bD]};        %Bounds for each parameters;
%
%% OUTPUTS:
% MLE:         maximum likelihood estimate
% S_laplace:   Laplace approximation of the covariance matrix evaluated at
%              MLE. Note: The Laplace approximation is computed in the
%              transformed space when bounds are provided
% x_NR:        Newton-Raphson intermediary states
% fx_NR:       function value corresponding to Newton-Raphson intermediary
%              states
%%
function [MLE,S_laplace,x_NR,fx_NR]=Newton_Raphson(fx,x0,varargin)
args = varargin;
nargs = length(args);
%% Default values
display='yes';              %Display process comments
output = 'transform';       %perform laplace approximation in the original space
laplace = 'yes';            %perform laplace approximation yes/no
log_transform = 'yes';      %perform Newton-Raphson in the log-space
nb_failed_iteration_limit=3;%maximal allowable number of failed iterations
convergence_tol=1E-4;       %convergence tolerance
delta = 1E-3;               %Relative delta x for numeric derivatives
bounds = [];                %Bounds for each parameters;

%% If provided, employ user-specific arguments
for i=1:2:nargs
    switch args{i}
        case 'output', output = args{i+1};
        case 'display', display = args{i+1};
        case 'laplace', laplace = args{i+1};
        case 'log_transform', log_transform = args{i+1};
        case 'convergence_tol', convergence_tol = args{i+1};
        case 'nb_failed_iteration_limit', nb_failed_iteration_limit = args{i+1};
        case 'delta', delta = args{i+1};
        case 'bounds', bounds = args{i+1};
        otherwise, error(['unrecognized argument ' args{i}])
    end
end

%% Analysis parameters
D=length(x0);               %number of parameters

if strcmp(display,'yes')
    disp(' ')
    disp('\Start Newton-Raphson algorithm')
    disp(['   \\Initial parameter values: ' num2str(x0)])
end
%% Transformation fct
if ~isempty(bounds); %apply transformation functions corresponding to the parameter space bounds
    [TR_fct,TR_fct_inv,dTR_fct,dTR_fct_inv]=transformation_functions(bounds);
    if strcmp(display,'yes')
        disp('   \\Transform parameter space ....')
    end
    if strcmp(log_transform,'yes')
        if strcmp(output,'transform')
            fx_TR=@(x_TR) log(abs(prod(dTR_fct_inv(x_TR)))*fx(TR_fct_inv(x_TR)));
        elseif strcmp(output,'original')
            fx_TR=@(x_TR) log(fx(TR_fct_inv(x_TR)));
        else error('Newton-Raphson: output is not valid')
        end
    else
        if strcmp(output,'transform')
            fx_TR=@(x_TR) sum(log(abs(dTR_fct_inv(x_TR))))+fx(TR_fct_inv(x_TR));
        elseif strcmp(output,'original')
            fx_TR=@(x_TR) fx(TR_fct_inv(x_TR));
        else error('Newton-Raphson: output is not valid')
        end
    end
    x0=TR_fct(x0);  %transform NR starting values into the transformed space
else %apply no transformation
    if strcmp(log_transform,'yes')
        fx_TR=@(x_TR) log(fx(x_TR));
    else
        fx_TR=@(x_TR) fx(x_TR);
    end
end

%% Gradient, 1D Hessian
%% double check for the objective computation 
%%
Term1=@(x,idx) fx_TR(x+idx*(max(delta*1E-4,delta*abs(x(logical(idx))))));
Term2=@(x,idx) fx_TR(x-idx*(max(delta*1E-4,delta*abs(x(logical(idx))))));
Term3=@(x,idx) max(delta*1E-4,delta*abs(x(logical(idx))));
grad_fct=@(x,idx) (Term1(x,idx)-Term2(x,idx))/(2*Term3(x,idx));                     %Gradient
hessian_1D_fct=@(x,idx) (Term1(x,idx)-2*fx_TR(x)+Term2(x,idx))/(Term3(x,idx)^2);    %Hessian
% grad_fct=@(x,idx) (fx_TR(x+idx*(max(delta*1E-4,delta*abs(x(logical(idx))))))-fx_TR(x-idx*max(delta*1E-4,delta*abs(x(logical(idx))))))/(2*max(delta*1E-4,delta*abs(x(logical(idx)))));                 %Gradient
% hessian_1D_fct=@(x,idx) (fx_TR(x+idx*max(delta*1E-4,delta*abs(x(logical(idx)))))-2*fx_TR(x)+fx_TR(x-idx*max(delta*1E-4,delta*abs(x(logical(idx))))))/((max(delta*1E-4,delta*abs(x(logical(idx)))))^2);%Hessian

%% Newton-Raphson algorithm
x_NR=x0;                    %Initialize parameter matrix
fx_NR=fx_TR(x0);            %Initialize fx vector
if strcmp(display,'yes')
    disp(['   \\Initial loglik: ' num2str(fx_NR)])
end

x_old=x_NR;                 %Initialize current states
fx_old=fx_TR(x_old);        %compute current function value
x_new=x_old;                %Initialize new proposed states
fx_new=[];                  %Initialize new state fx
hessian_NR=inf*zeros(1,D,'gpuArray');
loop=0;
param_fail_idx=zeros(1,D);
while any(hessian_NR>0)||isempty(fx_new)||(fx_NR(end)-fx_NR(end-1))>convergence_tol
    loop=loop+1;
    if strcmp(display,'yes')
        disp(['   \\loop #' num2str(loop)])
    end
    for d=1:D
        lambda=1;           %Initialize search step size
        idx=zeros(1,D);     %Initialize parameter index
        idx(d)=1;           %select current parameter index
        grad_old=grad_fct(x_old,idx);          %compute gradient
        hessian_old=hessian_1D_fct(x_old,idx); %compute 1D Hessian
        hessian_NR(d)=hessian_old;
        lambda_loop=0;      %initialize lambda loop count
        while any(isnan(x_new))||any(isinf(x_new))||isempty(fx_new)||fx_new<=fx_old %loop until the function value increase for the proposed state
            if lambda_loop>nb_failed_iteration_limit;
                break       %stop optimizing the current parameter if optimization fail > nb_failed_iteration_limit
            end
            if grad_old~=0&&hessian_old~=0&&~isnan(grad_old)&&~isnan(hessian_old)
                x_new(d)=x_old(d)-lambda*gather(grad_old/-abs(hessian_old));
            end
            try
                fx_new=fx_TR(x_new);
                
                %if isnan(fx_new)||isinf(fx_new)
                %    disp('c')
                %end
            catch err
                fx_new=-inf; %The likelihood cannot be evaluated for this parameter value
            end
            lambda=lambda/2;
            lambda_loop=lambda_loop+1;
        end
        if lambda_loop<=nb_failed_iteration_limit
            x_old=x_new;
            fx_old=fx_new;
            fx_NR=[fx_NR;fx_new];
            x_NR=[x_NR;x_new];
        else
            param_fail_idx(d)=1;
        end
        if all(param_fail_idx)
            break
        end
    end
    if strcmp(display,'yes')
        if strcmp(output,'original')
            MLE=TR_fct_inv(x_NR(end,:)); %Maximum likelihood estimate obtain with NR
        else
            MLE=x_NR(end,:); %Maximum likelihood estimate obtain with NR
        end
        disp(['     loglik=' num2str(fx_NR(end))])
        disp(['     parameters:' num2str(MLE) ])
    end
    if loop>100*nb_failed_iteration_limit;
        if strcmp(display,'yes')
            disp(['   \\\The number of failed iteration exceeds the allowable limit -> stop'])
        end
        break       %stop optimizing the current parameter if optimization fail > nb_failed_iteration_limit
    end
    if all(param_fail_idx)
        if strcmp(display,'yes')
            disp(['   \\\All parameters have failed to be optimized -> stop'])
        end
        break
    end
end
if strcmp(output,'original')
    MLE=TR_fct_inv(x_NR(end,:)); %Maximum likelihood estimate obtain with NR
else
    MLE=x_NR(end,:); %Maximum likelihood estimate obtain with NR
end
if strcmp(display,'yes')
    disp(' ')
    disp(['   \\Optimized parameters: ' num2str(MLE)])
    disp(' ')
end
%% Laplace approximation in the transformed space
if strcmp(laplace,'yes')&&strcmp(output,'transform')
    if strcmp(display,'yes')
        disp('   \\Compute Laplace approximation of the posterior in the transformed space')
    end
    H=numerical_hessian(MLE,fx_TR);%compute hessian matrix
elseif strcmp(laplace,'yes')&&strcmp(output,'original')
    if strcmp(display,'yes')
        disp('   \\Compute Laplace approximation of the posterior in the original space')
    end
    H=numerical_hessian(MLE,fx);%compute hessian matrix
else
    S_laplace=[];
end
if ~strcmp(laplace,'no')
    S_laplace=-(H^-1);                  %compute the laplace approximation of the posterior covariance
    S_laplace=(S_laplace+S_laplace')/2; %average the matrix and its transpose to make ensure numerical error does not affect the matrix symmetry
    try chol(S_laplace);
    catch err
        D=diag(sqrt(abs((diag(diag(H))^-1))));                  %compute the laplace approximation of the posterior covariance
        R=zeros(length(D));
        for i=1:length(D)
            for j=i+1:length(D)
                S_ij=-(H([i,j],[i,j])^-1);
                R(i,j)=sign(S_ij(1,1))*sign(S_ij(2,2))*S_ij(1,2)/sqrt(abs(S_ij(1,1)))/sqrt(abs(S_ij(2,2)));
            end
        end
        R=R+R'+eye(length(D));
        S_laplace=psd_constrain(diag(D)*R*diag(D));
        disp('        Warning: Newton_Raphson.m - Laplace approximation is not PSD')
        disp('        S_laplace = PSD_constrain(S_laplace)')
    end
end
%% Transform MLE back into the original space
if ~isempty(bounds);
    x_NR=TR_fct_inv(x_NR);
end