if IncludeStructuralAtt
    % TAGI parameters
    StructuralAttributes=0;
    TrueSpeedValues=[];
    if isempty(AnnModel)
        AnnModel.batch_size = 2;
        AnnModel.max_epoches = 30;
        AnnModel.cont_training = 0;
        AnnModel.evaluate_mode = 1;
        AnnModel.theta_nn = [];
    end
else
    AnnModel=[];
    RegressionModel=[];
    Stored_RegressionModel{model_i, elem_j}=[];
end

% Factoring Structural Attributes
% IncludeStructuralAtt

if ~isempty(AnnModel) && IncludeStructuralAtt
    StructuralAttributes=ElementData.StrucAtt;
    Stored_AnnModel = AnnModel;
    KRparam=zeros(2,1);
    KernelParameters{1}=0;
    KernelParameters{2}=0;
    InitialEx=[];
    InitialVar=[];
elseif ~IncludeStructuralAtt
    KRparam=zeros(2,1);
    KernelParameters{1}=0;
    KernelParameters{2}=0;
    InitialEx=[];
    InitialVar=[];
    Stored_RegressionModel{model_i, elem_j}=[];
    AnnModel=[];
end