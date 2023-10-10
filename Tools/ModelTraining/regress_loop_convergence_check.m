    if RegressLoop>1
        if SumLL(RegressLoop,1) > gather(sum(LoglikStoreValidation(:)))%&& ~stall_reached %SumLL(RegressLoop-1,1)
            AnnModel_Store=AnnModel;
            LoglikStoreValidation=LogLikValidation;
            RegressLoop_Store = RegressLoop;
            Best_RegressLoop = RegressLoop;
            fprintf('LL improved...\n');
        else
            Converge=1;
            AnnModel = AnnModel_Store;
            LogLikValidation=LoglikStoreValidation;
            Best_RegressLoop = RegressLoop_Store;
            fprintf('LL did not improve...\n');
        end
    % store the better of first iteration or the loaded model
    elseif SumLL(RegressLoop,1) > LogLikValidation_loaded
        disp('loaded AnnModel NOT stored...')
        AnnModel_Store=AnnModel;
        LoglikStoreValidation=LogLikValidation;
        RegressLoop_Store = RegressLoop;
        Best_RegressLoop = RegressLoop;
    else
        disp('loaded AnnModel stored... terminating regress loop')
        Converge=1;
        AnnModel=AnnModel_loaded;
        LogLikValidation=LogLikValidation_loaded;
        Best_RegressLoop = 0;
    end