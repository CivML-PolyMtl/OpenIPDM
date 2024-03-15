function [mu_a, variance_a, cov_ah] = mixture_truncation(x, p, cov, constrained_variable, l, u)
    % Perform full truncation and compute the mean and variance.

    mu_z = constrained_variable * x;
    sigma_z = sqrt(constrained_variable * p * constrained_variable');

    z1 = (-mu_z + u) / sigma_z;
    z2 = (mu_z - l) / sigma_z;

    Phi1 = normcdf(z1);
    Phi2 = normcdf(z2);

    phi1 = normpdf(z1);
    phi2 = normpdf(z2);

    mu_l = (mu_z - l) * Phi2 + sigma_z * phi2 + l;
    mu_u = (mu_z - u) * Phi1 - sigma_z * phi1 + u;
    mu_a = mu_l + mu_u - mu_z;

    term1 = ((-mu_z + u)^2 + sigma_z^2) * Phi1 + (-mu_z + u) * sigma_z * phi1;
    term2 = ((mu_z - l)^2 + sigma_z^2) * Phi2 + (mu_z - l) * sigma_z * phi2;
    term3 = -sigma_z^2 - (mu_z - mu_a)^2;
    term4 = (mu_z - ((mu_z - l) * Phi2 + sigma_z * phi2 + l))^2 + (mu_z - ((mu_z - u) * Phi1 - sigma_z * phi1 + u))^2;

    variance_a = term1 - (mu_u - u)^2 + term2 - (mu_l - l)^2 + term3 + term4;
    cov_az = term1 + (mu_u - u)*(-mu_z + u) + term2 - (mu_l - l)*(mu_z - l) - sigma_z^2;

    cov_ah = (cov_az/sigma_z^2) * cov;
    cov_ah = reshape(cov_ah, [], 1);
end

