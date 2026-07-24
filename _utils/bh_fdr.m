function p_fdr = bh_fdr(p)
% Benjamini-Hochberg FDR correction (no toolbox required)
    n = length(p);
    [p_sorted, sort_idx] = sort(p(:));
    p_fdr_sorted = p_sorted .* n ./ (1:n)';
    for k = n-1:-1:1
        p_fdr_sorted(k) = min(p_fdr_sorted(k), p_fdr_sorted(k+1));
    end
    p_fdr           = nan(size(p));
    p_fdr(sort_idx) = min(p_fdr_sorted, 1);
end