"""
stats_utils.py
--------------
Shared statistical helper functions for contingency analysis and
multiple-comparison correction.

Functions
---------
chi2_contingency_manual(cont)
    Chi-square test on a 2D contingency table (no external toolbox).

bh_fdr(p_vals)
    Benjamini–Hochberg false discovery rate (FDR) correction.

odds_ratio_woolf(a, b, c, d)
    Odds ratio and 95% confidence interval using Woolf's logit method.

Steven Errington, 2026
"""

from __future__ import annotations

import numpy as np
from scipy import stats


# ── Contingency table statistics ─────────────────────────────────────────────

def chi2_contingency_manual(cont: np.ndarray):
    """
    Perform a chi-square test on a 2D contingency table.

    This implementation computes expected counts and the chi-square
    statistic directly, without relying on ``scipy.stats.chi2_contingency``.

    Parameters
    ----------
    cont : np.ndarray
        2-D contingency table of observed counts.

    Returns
    -------
    chi2 : float
        Chi-square test statistic.
    p : float
        p-value from the chi-square distribution.
    df : int
        Degrees of freedom.

    Notes
    -----
    Expected counts are computed as the outer product of row and column
    sums divided by the total sample size.
    """
    cont  = cont.astype(float)
    row_s = cont.sum(axis=1, keepdims=True)
    col_s = cont.sum(axis=0, keepdims=True)
    N     = cont.sum()

    exp   = (row_s @ col_s) / N
    chi2  = np.nansum((cont - exp) ** 2 / exp)

    df    = (cont.shape[0] - 1) * (cont.shape[1] - 1)
    p     = 1 - stats.chi2.cdf(chi2, df)

    return chi2, p, df


# ── Multiple comparisons ─────────────────────────────────────────────────────

def bh_fdr(p_vals: np.ndarray) -> np.ndarray:
    """
    Apply Benjamini–Hochberg FDR correction.

    Controls the expected proportion of false discoveries among rejected
    hypotheses.

    Parameters
    ----------
    p_vals : np.ndarray
        Array of p-values.

    Returns
    -------
    np.ndarray
        FDR-adjusted p-values (same shape as input).

    Notes
    -----
    The procedure sorts p-values, applies the BH scaling factor, and
    enforces monotonicity via a cumulative minimum from the right.
    """
    p_vals = np.asarray(p_vals)
    n      = len(p_vals)

    sort_idx   = np.argsort(p_vals)
    p_sorted   = p_vals[sort_idx]
    p_fdr_sort = p_sorted * n / (np.arange(n) + 1)

    # Enforce monotonicity
    for k in range(n - 2, -1, -1):
        p_fdr_sort[k] = min(p_fdr_sort[k], p_fdr_sort[k + 1])

    p_fdr = np.empty(n)
    p_fdr[sort_idx] = np.minimum(p_fdr_sort, 1.0)

    return p_fdr


# ── Effect size ──────────────────────────────────────────────────────────────

def odds_ratio_woolf(a: float, b: float, c: float, d: float):
    """
    Compute odds ratio and 95% confidence interval (Woolf method).

    Uses the logit (log odds) approximation to estimate the confidence
    interval. Applies the Haldane–Anscombe correction if any cell count
    is zero.

    Parameters
    ----------
    a, b, c, d : float
        Cell counts of a 2x2 contingency table:

            [[a, b],
             [c, d]]

    Returns
    -------
    or_val : float
        Estimated odds ratio.
    ci_lo : float
        Lower bound of the 95% confidence interval.
    ci_hi : float
        Upper bound of the 95% confidence interval.

    Notes
    -----
    The standard error is computed on the log scale:
        SE = sqrt(1/a + 1/b + 1/c + 1/d)
    """
    if 0 in (a, b, c, d):
        a, b, c, d = a + 0.5, b + 0.5, c + 0.5, d + 0.5

    log_or    = np.log((a * d) / (b * c))
    se_log_or = np.sqrt(1/a + 1/b + 1/c + 1/d)

    or_val = np.exp(log_or)
    ci_lo  = np.exp(log_or - 1.96 * se_log_or)
    ci_hi  = np.exp(log_or + 1.96 * se_log_or)

    return or_val, ci_lo, ci_hi