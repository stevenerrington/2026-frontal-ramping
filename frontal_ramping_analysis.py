#!/usr/bin/env python3
"""Standalone frontal ramping analysis script.

This script reproduces the notebook workflow from ``frontal_ramping_analysis.ipynb``.
It loads the spike/session log, computes fixation-aligned SDFs by foreperiod
condition, z-scores the traces, detects ramping neurons, computes prevalence
statistics, and saves summary figures.
"""

from __future__ import annotations

import argparse
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy.stats as stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm.auto import tqdm

from _utils.analysis_utils import smooth, load_beh_mat, load_sdf_mat
from _utils.plot_utils import sig_marker, add_region_dividers
from _utils.stats_utils import bh_fdr, chi2_contingency_manual, odds_ratio_woolf


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description='Standalone frontal ramping analysis script.'
    )
    parser.add_argument(
        '--root', type=Path,
        default=Path(__file__).resolve().parent,
        help='Project root directory (default: script directory).'
    )
    parser.add_argument(
        '--raw-data', type=Path,
        default=Path('/Volumes/Mnemosyne/Data/2026_macaque_value/spk/'),
        help='Directory containing raw _spk.mat behavioural files.'
    )
    parser.add_argument(
        '--proc-data', type=Path,
        default=Path('/Volumes/Mnemosyne/Data/2026_macaque_value/proc/'),
        help='Directory containing processed SDF .mat files.'
    )
    parser.add_argument(
        '--recompute', action='store_true',
        help='Force recomputation of fixation arrays and linearity analysis.'
    )
    parser.add_argument(
        '--skip-figures', action='store_true',
        help='Skip figure generation and only compute data outputs.',
    )
    parser.add_argument(
        '--verbose', action='store_true',
        help='Print progress and diagnostic information.',
    )
    return parser.parse_args()


def get_dirs(root: Path, raw_data: Path, proc_data: Path) -> dict[str, Path]:
    dirs = {
        'root': root,
        'raw_data': raw_data,
        'proc_data': proc_data,
        'figures': root / '_figures',
    }
    data_dir = root / '_data'
    data_dir.mkdir(parents=True, exist_ok=True)
    dirs['figures'].mkdir(parents=True, exist_ok=True)
    dirs['data'] = data_dir
    return dirs


def compute_fixation_array(
    spike_log: pd.DataFrame,
    dirs: dict[str, Path],
    smooth_win: int = 100,
    force_recompute: bool = False,
    verbose: bool = False,
) -> np.ndarray:
    save_path = dirs['data'] / 'fixation_array.npy'
    if save_path.exists() and not force_recompute:
        if verbose:
            print('Loading existing fixation_array from disk...')
        return np.load(save_path)

    first_sdf_path = dirs['proc_data'] / f"{spike_log['neuron_label'].iloc[0]}.mat"
    first_sdf = load_sdf_mat(first_sdf_path)
    n_timepoints = first_sdf['sdf']['stim_on'].shape[1]

    n_neurons = len(spike_log)
    fixation_array = np.full((n_neurons, n_timepoints, 4), np.nan)

    if verbose:
        print(f'Computing fixation_array for {n_neurons} neurons...')

    foreperiod_masks = [
        ('all', lambda fp: fp > 0.5),
        ('short', lambda fp: (fp > 0.70) & (fp < 0.75)),
        ('medium', lambda fp: (fp > 0.90) & (fp < 1.10)),
        ('long', lambda fp: fp > 1.20),
    ]

    for neuron_i, row in tqdm(spike_log.iterrows(), total=n_neurons):
        beh_path = dirs['raw_data'] / f"{row['session']}_spk.mat"
        sdf_path = dirs['proc_data'] / f"{row['neuron_label']}.mat"
        try:
            beh = load_beh_mat(beh_path)
            sdf = load_sdf_mat(sdf_path)['sdf']['stim_on']
        except Exception as exc:
            print(f'  [WARN] neuron {neuron_i}: {exc}')
            continue

        foreperiod = beh['stim_on'] - beh['fixcross_fix']

        for cond_i, (_, cond_mask_fn) in enumerate(foreperiod_masks):
            mask = cond_mask_fn(foreperiod)
            if mask.sum() == 0:
                continue
            trial_mean = np.nanmean(sdf[mask, :], axis=0)
            fixation_array[neuron_i, :, cond_i] = smooth(trial_mean, smooth_win)

    np.save(save_path, fixation_array)
    if verbose:
        print(f'fixation_array saved to {save_path}')
    return fixation_array


def zscore_fixation_array(fixation_array: np.ndarray) -> np.ndarray:
    z_fixation_array = np.zeros_like(fixation_array)
    for neuron_i in range(fixation_array.shape[0]):
        for cond_i in range(fixation_array.shape[2]):
            trace = fixation_array[neuron_i, :, cond_i]
            mu = np.nanmean(trace)
            sd = np.nanstd(trace, ddof=1)
            if sd > 0:
                z_fixation_array[neuron_i, :, cond_i] = (trace - mu) / sd
    return z_fixation_array


def compute_rsq(y: np.ndarray, t: np.ndarray) -> float:
    coeffs = np.polyfit(t, y, 1)
    y_hat = np.polyval(coeffs, t)
    ss_res = np.sum((y - y_hat) ** 2)
    ss_tot = np.sum((y - y.mean()) ** 2)
    if ss_tot == 0:
        return 0.0
    return 1.0 - ss_res / ss_tot


def detect_ramping_neurons(
    z_fixation_array: np.ndarray,
    spike_log: pd.DataFrame,
    data_dir: Path,
    n_perm: int = 1000,
    p_thresh: float = 0.001,
    force_recompute: bool = False,
    verbose: bool = False,
) -> pd.DataFrame:
    linearity_path = data_dir / 'linearity_analysis.csv'
    if linearity_path.exists() and not force_recompute:
        if verbose:
            print('Loading existing linearity_analysis.csv...')
        return pd.read_csv(linearity_path)

    SDF_ZERO = 2000
    RAMP_MS = np.arange(-750, 1)
    RAMP_IDX = SDF_ZERO + RAMP_MS
    N_PTS = len(RAMP_MS)
    t = RAMP_MS.astype(float)

    rng = np.random.default_rng(42)
    records: list[dict[str, float | int]] = []
    n_neurons = z_fixation_array.shape[0]

    if verbose:
        print(f'Running ramp detection for {n_neurons} neurons ({n_perm} permutations each)...')

    for neuron_i in tqdm(range(n_neurons)):
        fr = z_fixation_array[neuron_i, RAMP_IDX, 0]
        if np.all(np.isnan(fr)):
            records.append({'flag': 0, 'pval': np.nan, 'slope': np.nan, 'r2': np.nan})
            continue

        slope = np.polyfit(t, fr, 1)[0]
        rsq = compute_rsq(fr, t)

        rsq_null = np.empty(n_perm)
        for perm_i in range(n_perm):
            shift_amt = rng.integers(1, N_PTS)
            fr_shift = np.roll(fr, shift_amt)
            rsq_null[perm_i] = compute_rsq(fr_shift, t)

        p_val = np.mean(rsq_null >= rsq)
        records.append({'flag': int(p_val < p_thresh), 'pval': p_val, 'slope': slope, 'r2': rsq})

    linearity_df = pd.DataFrame(records)
    linearity_df['area'] = spike_log['area'].values
    linearity_df['neuron_label'] = spike_log['neuron_label'].values
    linearity_df.to_csv(linearity_path, index=False)

    if verbose:
        print(f'linearity_analysis.csv saved to {linearity_path}')
    return linearity_df


def classify_ramping_neurons(linearity_df: pd.DataFrame, spike_log: pd.DataFrame) -> pd.DataFrame:
    R2_THRESH = 0.8
    pos_mask = (linearity_df['flag'] == 1) & (linearity_df['slope'] > 0) & (linearity_df['r2'] > R2_THRESH)
    neg_mask = (linearity_df['flag'] == 1) & (linearity_df['slope'] < 0) & (linearity_df['r2'] > R2_THRESH)

    analysis_log = spike_log.copy()
    analysis_log['ramping_fp_flag'] = 0
    analysis_log.loc[pos_mask, 'ramping_fp_flag'] = 1
    analysis_log.loc[neg_mask, 'ramping_fp_flag'] = -1
    return analysis_log


def compute_prevalence_and_stats(
    analysis_log: pd.DataFrame,
    linearity_df: pd.DataFrame,
    spike_log: pd.DataFrame,
    min_n: int = 10,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    areas = sorted(analysis_log['area'].unique())
    rows: list[dict[str, float | int | str]] = []
    for area in areas:
        flags = analysis_log.loc[analysis_log['area'] == area, 'ramping_fp_flag']
        n_total = len(flags)
        n_pos = int((flags == 1).sum())
        n_neg = int((flags == -1).sum())
        n_nonramp = int((flags == 0).sum())
        rows.append({
            'area': area,
            'n_total': n_total,
            'n_pos': n_pos,
            'n_neg': n_neg,
            'n_nonramp': n_nonramp,
            'n_ramp': n_pos + n_neg,
            'pct_ramp': 100 * (n_pos + n_neg) / n_total if n_total > 0 else np.nan,
        })

    prevalence_df = pd.DataFrame(rows).sort_values('pct_ramp', ascending=False)

    AREA_ORDER = [
        '24c',
        '6DR', '6DC', '6VaVb',
        '8B', '8A', '46d', '46df', '46v',
        '44', '45',
        '12r', '12m', '12o', '12l',
        'AI',
        '13l', '13m', '11ml',
        'cd', 'pu',
        'AMG',
    ]

    REGION_MAP_FULL = {
        '24c': ('MFC', '#993399'),
        '6DR': ('PMC', '#E67F14'), '6DC': ('PMC', '#E67F14'), '6VaVb': ('PMC', '#E67F14'),
        '8B': ('dlPFC', '#339EE5'), '8A': ('dlPFC', '#339EE5'), '46d': ('dlPFC', '#339EE5'),
        '46df': ('dlPFC', '#339EE5'), '46v': ('dlPFC', '#339EE5'),
        '44': ('IFG', '#4DB366'), '45': ('IFG', '#4DB366'),
        '12r': ('vlPFC', '#1A66B3'), '12m': ('vlPFC', '#1A66B3'),
        '12o': ('vlPFC', '#1A66B3'), '12l': ('vlPFC', '#1A66B3'),
        'AI': ('AI', '#B3B34C'),
        '13l': ('OFC', '#CC4C4C'), '13m': ('OFC', '#CC4C4C'), '11ml': ('OFC', '#CC4C4C'),
        'cd': ('Striatum', '#808080'), 'pu': ('Striatum', '#808080'),
        'AMG': ('AMG', '#4C4C4C'),
    }

    pt = prevalence_df[prevalence_df['n_total'] >= min_n].copy()
    order_present = [a for a in AREA_ORDER if a in pt['area'].str.strip().values]
    pt = pt.set_index('area').loc[order_present].reset_index()

    for df in [pt, prevalence_df, analysis_log, linearity_df, spike_log]:
        if 'area' in df.columns:
            df['area'] = df['area'].str.strip()

    pt['region'] = pt['area'].map(lambda a: REGION_MAP_FULL.get(a, ('Unknown', '#AAAAAA'))[0])
    pt['color'] = pt['area'].map(lambda a: REGION_MAP_FULL.get(a, ('Unknown', '#AAAAAA'))[1])

    overall_ramp = pt['n_ramp'].sum()
    overall_total = pt['n_total'].sum()
    overall_nonramp = overall_total - overall_ramp

    chi2_stat_list, chi2_p_list = [], []
    or_list, ci_lo_list, ci_hi_list = [], [], []

    for _, row in pt.iterrows():
        a = float(row['n_ramp'])
        b = float(row['n_nonramp'])
        c = float(overall_ramp - a)
        d = float(overall_nonramp - b)

        cont = np.array([[a, b], [c, d]], dtype=float)
        row_s = cont.sum(axis=1, keepdims=True)
        col_s = cont.sum(axis=0, keepdims=True)
        N_ = cont.sum()
        exp_ = (row_s @ col_s) / N_
        yates = np.sum((np.abs(cont - exp_) - 0.5) ** 2 / exp_)
        p_yates = 1 - stats.chi2.cdf(yates, 1)

        if (exp_ < 5).any():
            warnings.warn(f'Area {row["area"]}: expected cell < 5')

        chi2_stat_list.append(yates)
        chi2_p_list.append(p_yates)

        or_val, ci_lo, ci_hi = odds_ratio_woolf(a, b, c, d)
        or_list.append(or_val)
        ci_lo_list.append(ci_lo)
        ci_hi_list.append(ci_hi)

    chi2_p_arr = np.array(chi2_p_list)
    chi2_p_fdr = bh_fdr(chi2_p_arr)

    pt = pt.copy()
    pt['chi2_stat'] = chi2_stat_list
    pt['chi2_p_raw'] = chi2_p_arr
    pt['chi2_p_fdr'] = chi2_p_fdr
    pt['sig_fdr'] = chi2_p_fdr < 0.05
    pt['odds_ratio'] = or_list
    pt['or_ci_lo'] = ci_lo_list
    pt['or_ci_hi'] = ci_hi_list

    return pt, prevalence_df


def save_prevalence_figures(
    pt: pd.DataFrame,
    dirs: dict[str, Path],
    chi2_omni_ramp: float,
    p_omni_ramp: float,
    df_omni_ramp: int,
    chi2_omni_pn: float,
    p_omni_pn: float,
) -> None:
    pos_pct = 100 * pt['n_pos'] / pt['n_total']
    neg_pct = 100 * pt['n_neg'] / pt['n_total']
    x = np.arange(len(pt))
    colors = pt['color'].values

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1]}, constrained_layout=True)

    ax = axes[0]
    for xi, pos, neg, col in zip(x, pos_pct, neg_pct, colors):
        ax.bar(xi, pos, color=col, width=0.7, linewidth=0)
        ax.bar(xi, -neg, color=col, alpha=0.5, width=0.7, linewidth=0)

    ax.axhline(0, color='k', linewidth=0.8)
    add_region_dividers(ax, pt)
    ax.set_xlim(-0.5, len(pt) - 0.5)
    ax.set_ylim(-30, 30)
    ax.set_xticks([])
    ax.set_yticks(range(-30, 31, 10))
    ax.set_yticklabels([str(abs(v)) for v in range(-30, 31, 10)])
    ax.set_ylabel('Ramp prevalence (%)')
    ax.text(-1.2, 15, 'Positive', rotation=90, ha='center', va='center', fontsize=8)
    ax.text(-1.2, -15, 'Negative', rotation=90, ha='center', va='center', fontsize=8)

    title_str = (
        f'Omnibus ramp vs non-ramp: X²({df_omni_ramp})={chi2_omni_ramp:.1f}, p={p_omni_ramp:.4f}\n'
        f'Omnibus pos vs neg: X²({df_omni_pn})={chi2_omni_pn:.1f}, p={p_omni_pn:.4f}'
    )
    ax.set_title(title_str, fontsize=8, loc='left')
    sns.despine(ax=ax, bottom=True)

    ax2 = axes[1]
    ax2.axhline(1, color='k', linewidth=0.8, linestyle='--')
    add_region_dividers(ax2, pt)

    for i, row in pt.reset_index(drop=True).iterrows():
        col = row['color']
        or_ = row['odds_ratio']
        ci_lo = row['or_ci_lo']
        ci_hi = row['or_ci_hi']
        ax2.plot([i, i], [ci_lo, ci_hi], color=col, linewidth=1.5)
        if row['sig_fdr']:
            ax2.scatter(i, or_, s=60, color=col, zorder=5)
        else:
            ax2.scatter(i, or_, s=60, facecolors='none', edgecolors=col, linewidths=1.2, zorder=5)
        marker = sig_marker(row['chi2_p_fdr'])
        if marker:
            ax2.text(i, 12, marker, ha='center', va='center', fontsize=8)
        ax2.text(i, 0.12, f"n={int(row['n_total'])}", ha='center', va='center', fontsize=6)

    ax2.set_yscale('log')
    ax2.set_ylim(0.05, 20)
    ax2.set_xlim(-0.5, len(pt) - 0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(pt['area'], rotation=45, ha='right', fontsize=9)
    ax2.set_ylabel('Odds ratio vs rest\n(log scale)')
    sns.despine(ax=ax2)

    fname = dirs['figures'] / 'foreperiod_prevalence_and_OR.pdf'
    fig.savefig(fname, bbox_inches='tight')
    plt.close(fig)


def save_foreperiod_sdf_figures(
    analysis_log: pd.DataFrame,
    z_fixation_array: np.ndarray,
    pt: pd.DataFrame,
    dirs: dict[str, Path],
) -> None:
    TIME_VEC = np.arange(-2000, 2001)
    PLOT_MASK = (TIME_VEC >= -1500) & (TIME_VEC <= 250)
    TIME_PLOT = TIME_VEC[PLOT_MASK]
    FP_TIMES = [-1300, -1000, -700]

    COND_LABELS = ['Short', 'Medium', 'Long']
    COND_COLORS = ['#4CA5E8', '#1A66B3', '#00173F']
    N_COLS = 4

    for flag, label in [(1, 'Positive ramp'), (-1, 'Negative ramp')]:
        areas_plot = pt['area'].values
        n_areas = len(areas_plot)
        n_rows = int(np.ceil(n_areas / N_COLS))

        fig, axes = plt.subplots(n_rows, N_COLS, figsize=(N_COLS * 3.5, n_rows * 2.8), constrained_layout=True)
        axes_flat = axes.flatten()
        fig.suptitle(f'Foreperiod SDF — {label} neurons', fontsize=13)

        for area_i, area in enumerate(areas_plot):
            ax = axes_flat[area_i]
            area_mask = (analysis_log['area'] == area) & (analysis_log['ramping_fp_flag'] == flag)
            neuron_indices = np.where(area_mask)[0]
            n_area = len(neuron_indices)

            ax.set_title(f'{area}  (n={n_area})', fontsize=9)
            ax.set_xlim(-1500, 250)
            ax.set_xlabel('Time from stim on (ms)', fontsize=7)
            ax.set_ylabel('FR (z)', fontsize=7)
            ax.tick_params(labelsize=7)
            sns.despine(ax=ax)

            if n_area == 0:
                ax.text(np.mean([-1500, 250]), 0, 'No neurons', ha='center', fontsize=8, color='#AAAAAA')
                continue

            area_data = z_fixation_array[neuron_indices, :, :][:, PLOT_MASK, :]
            for cond_i, (clabel, ccolor) in enumerate(zip(COND_LABELS, COND_COLORS)):
                cond_data = area_data[:, :, cond_i + 1]
                mu = np.nanmean(cond_data, axis=0)
                sem = np.nanstd(cond_data, axis=0, ddof=1) / np.sqrt(n_area)
                ax.fill_between(TIME_PLOT, mu - sem, mu + sem, color=ccolor, alpha=0.15)
                ax.plot(TIME_PLOT, mu, color=ccolor, linewidth=1.5, label=clabel)

            for fp in FP_TIMES:
                ax.axvline(fp, color='#BBBBBB', linewidth=0.8, linestyle='--')
            ax.axvline(0, color='#333333', linewidth=1.0)
            if area_i == 0:
                ax.legend(fontsize=7, loc='upper left', frameon=False)

        for ax in axes_flat[n_areas:]:
            ax.set_visible(False)

        fname = dirs['figures'] / f"foreperiod_sdf_{'pos' if flag == 1 else 'neg'}_ramp.pdf"
        fig.savefig(fname, bbox_inches='tight')
        plt.close(fig)


def build_feature_table(
    z_fixation_array: np.ndarray,
    analysis_log: pd.DataFrame,
) -> pd.DataFrame:
    SDF_ZERO = 2000
    WIN_MS = np.arange(-600, 1)
    WIN_IDX = SDF_ZERO + WIN_MS
    t_win = WIN_MS.astype(float)

    COND_NAMES = ['short', 'medium', 'long']
    COND_IDX = [1, 2, 3]

    records: list[dict[str, float | int | str]] = []
    ramp_mask = analysis_log['ramping_fp_flag'] != 0
    ramp_indices = np.where(ramp_mask)[0]

    for neuron_i in ramp_indices:
        row = analysis_log.iloc[neuron_i]
        ramp_flag = row['ramping_fp_flag']
        ramp_type = 'Positive' if ramp_flag == 1 else 'Negative'

        for cond_name, cond_i in zip(COND_NAMES, COND_IDX):
            trace = z_fixation_array[neuron_i, WIN_IDX, cond_i]
            if np.all(np.isnan(trace)):
                slope = np.nan
                mean_fr = np.nan
            else:
                slope, mean_fr = np.polyfit(t_win, trace, 1)[0], np.nanmean(trace)
            records.append({
                'neuron_i': neuron_i,
                'neuron_label': row['neuron_label'],
                'area': row['area'],
                'ramp_type': ramp_type,
                'condition': cond_name,
                'slope': slope,
                'mean_fr': mean_fr,
            })

    features_df = pd.DataFrame(records)
    features_df['condition'] = pd.Categorical(features_df['condition'], categories=COND_NAMES, ordered=True)
    return features_df


def main() -> None:
    args = parse_args()
    sns.set_theme(style='ticks', font_scale=0.9)
    plt.rcParams.update({
        'axes.spines.top': False,
        'axes.spines.right': False,
        'figure.dpi': 150,
    })

    dirs = get_dirs(args.root, args.raw_data, args.proc_data)
    spike_log = pd.read_csv(dirs['data'] / 'spike_log.csv')
    n_neurons = len(spike_log)

    if args.verbose:
        print(f'Loaded spike_log: {n_neurons} neurons across {spike_log["area"].nunique()} areas')

    fixation_array = compute_fixation_array(
        spike_log, dirs, smooth_win=100, force_recompute=args.recompute, verbose=args.verbose
    )
    if args.verbose:
        print(f'fixation_array shape: {fixation_array.shape}  (neurons x time x conditions)')

    z_fixation_array = zscore_fixation_array(fixation_array)
    if args.verbose:
        print('Z-scoring complete.')

    linearity_df = detect_ramping_neurons(
        z_fixation_array,
        spike_log,
        dirs['data'],
        n_perm=1000,
        force_recompute=args.recompute,
        verbose=args.verbose,
    )

    analysis_log = classify_ramping_neurons(linearity_df, spike_log)
    pt, prevalence_df = compute_prevalence_and_stats(analysis_log, linearity_df, spike_log)

    cont_ramp = pt[['n_ramp', 'n_nonramp']].values
    chi2_omni_ramp, p_omni_ramp, df_omni_ramp = chi2_contingency_manual(cont_ramp)

    ramp_rows = pt['n_ramp'] > 0
    cont_posneg = pt.loc[ramp_rows, ['n_pos', 'n_neg']].values
    chi2_omni_pn, p_omni_pn, df_omni_pn = chi2_contingency_manual(cont_posneg)

    if args.verbose:
        print(f'Omnibus ramp vs non-ramp: X²({df_omni_ramp}) = {chi2_omni_ramp:.2f}, p = {p_omni_ramp:.4f}')
        print(f'Omnibus pos vs neg:       X²({df_omni_pn}) = {chi2_omni_pn:.2f}, p = {p_omni_pn:.4f}')

    pt_path = dirs['data'] / 'prevalence_table.csv'
    pt.to_csv(pt_path, index=False)
    if args.verbose:
        print(f'Saved prevalence table to {pt_path}')

    if not args.skip_figures:
        save_prevalence_figures(
            pt, dirs, chi2_omni_ramp, p_omni_ramp, df_omni_ramp, chi2_omni_pn, p_omni_pn
        )
        save_foreperiod_sdf_figures(analysis_log, z_fixation_array, pt, dirs)
        if args.verbose:
            print(f'Saved figures to {dirs["figures"]}')

    features_df = build_feature_table(z_fixation_array, analysis_log)
    features_path = dirs['data'] / 'foreground_features.csv'
    features_df.to_csv(features_path, index=False)
    if args.verbose:
        print(f'Saved ramping feature table to {features_path}')


if __name__ == '__main__':
    main()
