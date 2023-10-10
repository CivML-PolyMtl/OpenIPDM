function [x_train, y_train, S_y_train, x_valid, y_valid, S_y_valid] = TAGI_train_test_split(X, y, Sy, split_ratio, seed, shuffle)
    if seed    
        rng(seed);
    else
        rng('shuffle');
    end

    if shuffle
        idxtrain  = randperm(size(y, 1));
        y         = y(idxtrain, :);
        Sy        = Sy(idxtrain, :);
        X         = X(idxtrain, :);
    end

    idxtrain = round(size(X, 1) * split_ratio);

    x_train   = X(1 : idxtrain, :);
    y_train   = y(1 : idxtrain, :);
    S_y_train = Sy(1 : idxtrain, :);
    x_valid   = X(idxtrain + 1 : end, :);  
    y_valid   = y(idxtrain + 1 : end, :);
    S_y_valid = Sy(idxtrain + 1 : end, :);
end