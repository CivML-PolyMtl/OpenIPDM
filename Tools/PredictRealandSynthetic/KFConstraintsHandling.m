function [xTrunc,PTrunc,Success]=KFConstraintsHandling(xTrunc,PTrunc,D,d,param)
r=length(d);Success=1;
for k = 1 :2: r-1
    try
        [Utrunc, Wtrunc, Vtrunc] = svd(PTrunc);
    catch
        sprintf('Variance has NAN values')
        Success=0;
        return;
    end
    Ttrunc = Utrunc;
    TTT = Ttrunc * Ttrunc';
    if (norm(eye(size(TTT)) - TTT) > 1e-8)
       disp('Error - Ttrunc is not orthogonal.');
       Success=0;
       return;
    end
    if (norm(Utrunc*Wtrunc*Utrunc' - PTrunc) > 1e-8)
        disp('Error - SVD failed for pdf trunction');
        Success=0;
        return;
    end
    % Compute the modified Gram-Schmidt transformation S * Amgs = [ Wmgs ; 0 ].
    % Amgs is a given n x m matrix, and S is an orthogonal n x n matrix, and Wmgs is an m x m matrix.
    Amgs = sqrt(Wtrunc) * Ttrunc' * D(k,:)'; % n x 1, where n = number of states
    [Wmgs, S] = MGS(Amgs);
    S = S * sqrt(D(k,:) * PTrunc * D(k,:)') / Wmgs;
    cTrunc = (d(k) - D(k,:) * xTrunc) / sqrt(D(k,:) * PTrunc * D(k,:)');
    dTrunc = (d(k+1) - D(k+1,:) * xTrunc) / sqrt(D(k+1,:) * PTrunc * D(k+1,:)');
    % The next 3 lines are for inequality constraints. In our example, they
    % are commented out because our problem uses equality constraints.
    alpha = sqrt(2/pi) / (erf(dTrunc/sqrt(2)) - erf(cTrunc/sqrt(2)));
    mu = alpha * (exp(-cTrunc^2/2) - exp(-dTrunc^2/2));
    sigma2 = alpha * (exp(-cTrunc^2/2) * (cTrunc - 2 * mu) - exp(-dTrunc^2/2) * (dTrunc - 2 * mu)) + mu^2 + 1;

    % The following two lines are used for equality constraints.
    % Since we have equality constraints, the mean of zTrunc = cTrunc = dTrunc,
    % and its variance is 0.
%             mu = cTrunc;
%             sigma2 = 0;
    zTrunc = zeros(size(xTrunc));
    CovZ = eye(length(zTrunc));
    if ~isnan(mu) %&& ~isinf(mu)
            zTrunc(1) = mu;
            CovZ(1,1) = sigma2;
    else
        sprintf('Constraints failed x=%d,xbar=%d',xTrunc(1),xTrunc(2))
        Success=0;
        return;
    end
    xTrunc = Ttrunc * sqrt(Wtrunc) * S' * zTrunc + xTrunc;
    PTrunc = Ttrunc * sqrt(Wtrunc) * S' * CovZ * S * sqrt(Wtrunc) * Ttrunc';
end
end