function [chi2_stat, p_val, df] = chi2_contingency(cont)
% Chi-square test on an arbitrary contingency matrix (no toolbox required)
    rowSums  = sum(cont, 2);
    colSums  = sum(cont, 1);
    N        = sum(cont(:));
    expected = (rowSums * colSums) / N;
    chi2_stat = sum(sum((cont - expected).^2 ./ expected));
    df        = (size(cont, 1) - 1) * (size(cont, 2) - 1);
    p_val     = 1 - chi2cdf(chi2_stat, df);
end
