"""
analysis_utils.py
----------------
Shared I/O, signal-processing, and neuron-metric utilities for the
frontal ramping analysis.

Functions
---------
smooth(x, window)
    Moving-average smoothing equivalent to MATLAB smooth(x, N).

load_sdf_mat(path)
    Load a processed SDF .mat file (v7 format).

load_beh_mat(path)
    Load a raw behavioural _spk.mat file (v7 format).

load_mat_hdf5(path)
    Load a processed proc .mat file (v7.3 / HDF5 format) containing
    both sdf and raster structs.

load_beh_hdf5(path)
    Load trial event times from a v7.3 _spk.mat file.

get_fp_masks(foreperiod)
    Return boolean trial masks for each foreperiod condition.

exp_decay(lag, A, tau_c, B)
    Exponential decay model for autocorrelation fitting.

compute_ff_cv(spike_times_list, win, scale)
    Compute Fano factor and CV from trial-level spike counts.

compute_timescale(sdf_fixation)
    Fit intrinsic timescale from fixation SDF autocorrelation.

process_neuron(args)
    Top-level worker function for parallel neuron metric computation.

Steven Errington, 2026
"""

from __future__ import annotations

import warnings
from pathlib import Path

import h5py
import numpy as np
import scipy.io as sio
from scipy.optimize import curve_fit


# ── Constants ─────────────────────────────────────────────────────────────────

FF_WIN_MS      = (-600, 0)       # pre-stimulus spike count window (ms)
FIX_BASELINE   = slice(0, 1000)  # -2000 to -1000 ms in fixation SDF
AC_LAG_MAX_MS  = 800             # maximum autocorrelation lag (ms)
MIN_TRIALS     = 5               # minimum trials for reliable FF/CV
MIN_MEAN_COUNT = 0.5             # minimum mean spike count per trial
TAU_BOUNDS     = (5, 2000)       # plausible timescale range (ms)
COND_NAMES     = ['all', 'short', 'medium', 'long']


# ── Signal processing ─────────────────────────────────────────────────────────

def smooth(x: np.ndarray, window: int) -> np.ndarray:
    """
    Moving-average smooth equivalent to MATLAB smooth(x, window).

    Uses a centred uniform kernel with shrinking windows at the edges,
    matching MATLAB's default 'moving' method behaviour.

    Parameters
    ----------
    x : np.ndarray
        1-D input signal.
    window : int
        Width of the smoothing window in samples.

    Returns
    -------
    np.ndarray
        Smoothed signal, same length as x.
    """
    if window < 2:
        return x.copy()
    kernel = np.ones(window) / window
    pad    = window // 2
    x_pad  = np.pad(x.astype(float), pad, mode='edge')
    out    = np.convolve(x_pad, kernel, mode='valid')
    return out[:len(x)]


# ── v7 I/O ────────────────────────────────────────────────────────────────────

def load_sdf_mat(path: Path | str) -> dict:
    """
    Load a processed SDF .mat file saved by get_spk_sdfs.m (v7 format).

    Parameters
    ----------
    path : Path or str
        Full path to the .mat file.

    Returns
    -------
    dict with key 'sdf', itself a dict of np.ndarray each shape
    (n_trials, n_timepoints): 'stim_on', 'fixation', 'reward', 'brk_fix'.

    Raises
    ------
    NotImplementedError
        If the file is v7.3 (HDF5). Use load_mat_hdf5 instead.
    """
    mat        = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    sdf_struct = mat['sdf']
    return {
        'sdf': {
            'stim_on':  np.atleast_2d(sdf_struct.stim_on),
            'fixation': np.atleast_2d(sdf_struct.fixation),
            'reward':   np.atleast_2d(sdf_struct.reward),
            'brk_fix':  np.atleast_2d(sdf_struct.brk_fix),
        },
    }


def load_beh_mat(path: Path | str) -> dict:
    """
    Load trial event times from a v7 _spk.mat file.

    Parameters
    ----------
    path : Path or str
        Full path to the _spk.mat file.

    Returns
    -------
    dict with keys 'stim_on', 'fixcross_fix', 'rew_on', 'brk_fix',
    each a np.ndarray of shape (n_trials,) in seconds.
    """
    mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    t   = mat['t_evt']
    return {
        'stim_on':      np.atleast_1d(t.stim_on),
        'fixcross_fix': np.atleast_1d(t.fixcross_fix),
        'rew_on':       np.atleast_1d(t.rew_on),
        'brk_fix':      np.atleast_1d(t.brk_fix),
    }


# ── v7.3 / HDF5 I/O ──────────────────────────────────────────────────────────

def load_mat_hdf5(path: Path | str) -> dict:
    """
    Load a v7.3 proc .mat file (HDF5) containing sdf and raster structs.

    MATLAB stores arrays column-major; h5py reads row-major, so SDF
    arrays are transposed to restore (n_trials, n_timepoints) layout.

    MATLAB cell arrays of spike times become HDF5 object reference arrays.
    Each reference is dereferenced to recover the per-trial spike vector.

    Parameters
    ----------
    path : Path or str
        Full path to the proc .mat file.

    Returns
    -------
    dict with keys:
        'sdf'    : dict of np.ndarray, shape (n_trials, n_timepoints)
        'raster' : dict of lists of np.ndarray (one array per trial)
    """
    with h5py.File(str(path), 'r') as f:
        sdf = {
            field: f['sdf'][field][:].T.astype(float)
            for field in ('stim_on', 'fixation', 'reward', 'brk_fix')
            if field in f['sdf']
        }
        raster = {}
        for event in ('stim_on', 'fixation', 'reward', 'brk_fix'):
            if event not in f['raster']:
                continue
            refs   = f['raster'][event][:].flatten()
            trials = []
            for ref in refs:
                spk = f[ref][()].flatten().astype(float)
                trials.append(spk if spk.size > 0 else np.array([]))
            raster[event] = trials

    return {'sdf': sdf, 'raster': raster}


def load_beh_hdf5(path: Path | str) -> np.ndarray:
    """
    Load foreperiod durations from a v7.3 _spk.mat file.

    Parameters
    ----------
    path : Path or str
        Full path to the _spk.mat file.

    Returns
    -------
    np.ndarray, shape (n_trials,)
        Foreperiod in seconds: stim_on - fixcross_fix.
    """
    with h5py.File(str(path), 'r') as f:
        stim_on      = f['t_evt']['stim_on'][:].flatten()
        fixcross_fix = f['t_evt']['fixcross_fix'][:].flatten()
    return stim_on - fixcross_fix


# ── Foreperiod masks ──────────────────────────────────────────────────────────

def get_fp_masks(foreperiod: np.ndarray) -> dict:
    """
    Return boolean trial masks for each foreperiod condition.

    Parameters
    ----------
    foreperiod : np.ndarray
        Per-trial foreperiod durations in seconds.

    Returns
    -------
    dict with keys 'all', 'short', 'medium', 'long', each a boolean
    np.ndarray of shape (n_trials,).
    """
    return {
        'short':  (foreperiod > 0.68) & (foreperiod < 0.75),
        'medium': (foreperiod > 0.95) & (foreperiod < 1.05),
        'long':   (foreperiod > 1.25) & (foreperiod < 1.35),
        'all': (
            ((foreperiod > 0.68) & (foreperiod < 0.75)) |
            ((foreperiod > 0.95) & (foreperiod < 1.05)) |
            ((foreperiod > 1.25) & (foreperiod < 1.35))
        )
    }
# ── Neuron metrics ────────────────────────────────────────────────────────────

def exp_decay(lag: np.ndarray, A: float, tau_c: float, B: float) -> np.ndarray:
    """
    Exponential decay model for autocorrelation fitting.

        AC(lag) = A * exp(-lag / tau_c) + B

    Parameters
    ----------
    lag : np.ndarray
        Lag values in ms.
    A : float
        Amplitude (bounded 0 to 2 in fitting).
    tau_c : float
        Timescale in ms.
    B : float
        Offset (bounded -1 to 1 in fitting).
    """
    return A * np.exp(-lag / tau_c) + B


def compute_ff_cv(
    spike_times_list: list,
    win: tuple = FF_WIN_MS,
    scale: float = 1.0,
) -> tuple[float, float]:
    """
    Compute Fano factor and CV from trial-level spike counts.

    Parameters
    ----------
    spike_times_list : list of np.ndarray
        Spike times per trial.
    win : tuple (t_start, t_end)
        Analysis window in ms.
    scale : float
        Multiply spike times by this to convert to ms.
        Use 1.0 if already ms, 1000.0 if seconds.

    Returns
    -------
    ff, cv : float
        Fano factor and coefficient of variation.
        Returns (nan, nan) if mean count < MIN_MEAN_COUNT or
        fewer than MIN_TRIALS trials.
    """
    t0, t1 = win
    counts = np.array([
        np.sum((spk * scale >= t0) & (spk * scale <= t1))
        for spk in spike_times_list
    ], dtype=float)

    if len(counts) < MIN_TRIALS or counts.mean() < MIN_MEAN_COUNT:
        return np.nan, np.nan

    ff = counts.var(ddof=1) / counts.mean()
    cv = counts.std(ddof=1) / counts.mean()
    return ff, cv


def compute_timescale(
    sdf_fixation: np.ndarray,
) -> tuple[float, float, float]:
    """
    Fit intrinsic timescale from fixation-aligned SDF autocorrelation.

    Uses the pre-task baseline window (-2000 to -1000 ms, FIX_BASELINE
    slice) and fits AC(lag) = A * exp(-lag / tau_c) + B over lags
    0 to AC_LAG_MAX_MS ms (Wasmuht et al. / Zeisler et al. approach).

    Parameters
    ----------
    sdf_fixation : np.ndarray, shape (n_trials, n_timepoints)
        SDF aligned to fixation cross, spanning -2000 to 2000 ms.

    Returns
    -------
    tau_c, A, B : float
        Fitted parameters. All NaN if fit fails or is out of bounds.
    """
    baseline = np.nanmean(sdf_fixation[:, FIX_BASELINE], axis=0)

    if np.all(np.isnan(baseline)) or baseline.std() < 1e-6:
        return np.nan, np.nan, np.nan

    b       = baseline - baseline.mean()
    ac_full = np.correlate(b, b, mode='full')
    ac_full = ac_full / ac_full[len(b) - 1]

    ac_pos   = ac_full[len(b) - 1:]
    lags     = np.arange(len(b), dtype=float)
    fit_mask = lags <= AC_LAG_MAX_MS

    try:
        popt, _ = curve_fit(
            exp_decay,
            lags[fit_mask], ac_pos[fit_mask],
            p0=(0.5, 100.0, 0.0),
            bounds=([0, TAU_BOUNDS[0], -1.0],
                    [2, TAU_BOUNDS[1],  1.0]),
            maxfev=5000,
        )
        return popt[1], popt[0], popt[2]   # tau_c, A, B

    except (RuntimeError, ValueError):
        return np.nan, np.nan, np.nan

def load_beh_mat_v7(path: Path | str) -> np.ndarray:
    """
    Load foreperiod from a standard v7 _spk.mat file via scipy.
    Returns stim_on - fixcross_fix in seconds.
    """
    mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    t   = mat['t_evt']
    return np.atleast_1d(t.stim_on) - np.atleast_1d(t.fixcross_fix)

# ── Parallel worker ───────────────────────────────────────────────────────────

def process_neuron(args: tuple) -> dict:
    """
    Compute all metrics for a single neuron.

    Designed as a top-level importable function so it can be pickled
    by ProcessPoolExecutor. Accepts a single tuple argument to allow
    use with executor.map().

    Parameters
    ----------
    args : tuple
        (neuron_i, neuron_label, area, session,
         sdf_path, beh_path, raster_scale)

    Returns
    -------
    dict
        Per-neuron metrics with keys:
            neuron_i, neuron_label, area, session,
            ff_{cond}, cv_{cond} for each cond in COND_NAMES,
            tau_c, tau_A, tau_B.
    """
    (neuron_i, neuron_label, area, session,
     sdf_path, beh_path, raster_scale) = args

    base = {
        'neuron_i':     neuron_i,
        'neuron_label': neuron_label,
        'area':         area,
        'session':      session,
        'tau_c': np.nan, 'tau_A': np.nan, 'tau_B': np.nan,
    }
    for cond in COND_NAMES:
        base[f'ff_{cond}'] = np.nan
        base[f'cv_{cond}'] = np.nan

    try:
        data = load_mat_hdf5(sdf_path)

        # ── Foreperiod masks ─────────────────────────────────────────────
        masks = None
        try:
            foreperiod = load_beh_mat_v7(beh_path)
            masks      = get_fp_masks(foreperiod)
        except Exception as e:
            warnings.warn(f'Beh load failed [{session}]: {e}')

        # ── FF and CV per foreperiod condition ───────────────────────────
        if masks is not None:
            spike_list = data['raster'].get('stim_on', [])
            n_raster   = len(spike_list)

            for cond, mask in masks.items():
                try:
                    trial_idx = np.where(mask)[0]
                    trial_idx = trial_idx[trial_idx < n_raster]
                    if len(trial_idx) == 0:
                        continue
                    subset = [spike_list[i] for i in trial_idx]
                    ff, cv = compute_ff_cv(subset, FF_WIN_MS, raster_scale)
                    base[f'ff_{cond}'] = ff
                    base[f'cv_{cond}'] = cv
                except Exception as e:
                    warnings.warn(
                        f'FF/CV ({cond}) failed [{neuron_label}]: {e}'
                    )

        # ── Intrinsic timescale ──────────────────────────────────────────
        try:
            sdf_fix       = np.atleast_2d(data['sdf']['fixation'])
            tau_c, A, B   = compute_timescale(sdf_fix)
            base['tau_c'] = tau_c
            base['tau_A'], base['tau_B'] = A, B
        except Exception as e:
            warnings.warn(f'Timescale failed [{neuron_label}]: {e}')

    except Exception as e:
        warnings.warn(f'Could not load [{neuron_label}]: {e}')

    return base
