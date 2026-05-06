"""
plot_utils.py
-------------
Shared plotting helper functions for statistical annotation and
visual segmentation of ordered data.

Functions
---------
sig_marker(p)
    Convert a p-value to significance asterisk notation.

add_region_dividers(ax, pt)
    Add dashed vertical lines at region boundaries on a plot.

Steven Errington, 2026
"""

from __future__ import annotations

import numpy as np
import pandas as pd


# ── Significance annotation ──────────────────────────────────────────────────

def sig_marker(p: float) -> str:
    """
    Convert a p-value to significance marker notation.

    Parameters
    ----------
    p : float
        p-value from a statistical test.

    Returns
    -------
    str
        Significance marker:
            '***' for p < 0.001
            '**'  for p < 0.01
            '*'   for p < 0.05
            ''    otherwise
    """
    if p < 0.001:
        return '***'
    if p < 0.01:
        return '**'
    if p < 0.05:
        return '*'
    return ''


# ── Plot annotation ──────────────────────────────────────────────────────────

def add_region_dividers(ax, pt: pd.DataFrame):
    """
    Add dashed vertical lines at boundaries between regions.

    Useful for plots where rows (e.g. neurons or channels) are sorted
    by cytoarchitectonic or anatomical region.

    Parameters
    ----------
    ax : matplotlib.axes.Axes
        Axis object to draw on.
    pt : pandas.DataFrame
        Table containing a ``'region'`` column specifying region labels
        in plotting order.

    Returns
    -------
    None

    Notes
    -----
    A vertical dashed line is drawn at each index where the region label
    changes between consecutive rows.
    """
    regions = pt['region'].values
    bounds  = np.where(regions[:-1] != regions[1:])[0]

    for b in bounds:
        ax.axvline(b + 0.5, color='#BBBBBB', linewidth=0.8, linestyle='--')