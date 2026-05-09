"""
ramping_utils.py
----------------
Shared I/O and signal-processing utilities for the frontal ramping analysis.

Functions
---------
smooth(x, window)
    Moving-average smoothing equivalent to MATLAB smooth(x, N).

load_sdf_mat(path)
    Load a processed SDF .mat file (output of get_spk_sdfs.m).

load_beh_mat(path)
    Load a raw behavioural _spk.mat file and return event times.

Steven Errington, 2026
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import scipy.io as sio
from scipy.optimize import curve_fit
import h5py


# ── Signal processing ────────────────────────────────────────────────────────

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


# ── I/O ─────────────────────────────────────────────────────────────────────

def load_sdf_mat(path: Path | str) -> dict:
    """
    Load a processed SDF .mat file saved by get_spk_sdfs.m.

    The file is expected to contain two MATLAB structs, ``sdf`` and
    ``raster``, each with fields ``stim_on``, ``fixation``, ``reward``,
    and ``brk_fix``.

    Parameters
    ----------
    path : Path or str
        Full path to the .mat file (e.g. ``proc/session_ch_clust.mat``).

    Returns
    -------
    dict with key ``'sdf'``, itself a dict:
        ``'stim_on'``  : np.ndarray, shape (n_trials, n_timepoints)
        ``'fixation'`` : np.ndarray, shape (n_trials, n_timepoints)
        ``'reward'``   : np.ndarray, shape (n_trials, n_timepoints)
        ``'brk_fix'``  : np.ndarray, shape (n_trials, n_timepoints)

    Raises
    ------
    NotImplementedError
        If the file was saved with MATLAB's ``-v7.3`` flag (HDF5 format).
        In that case use ``load_sdf_mat_hdf5`` instead.

    Notes
    -----
    ``scipy.io.loadmat`` is used with ``squeeze_me=True`` and
    ``struct_as_record=False`` so MATLAB struct fields are accessible
    as Python attributes (e.g. ``mat_struct.stim_on``).
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


def load_sdf_mat_hdf5(path: Path | str) -> dict:
    """
    Load a processed SDF .mat file saved with MATLAB's ``-v7.3`` flag.

    Use this as a drop-in replacement for ``load_sdf_mat`` if you get a
    ``NotImplementedError`` from scipy (v7.3 files use HDF5 internally).
    Requires ``h5py`` (``pip install h5py``).

    Parameters
    ----------
    path : Path or str
        Full path to the .mat file.

    Returns
    -------
    dict
        Same structure as ``load_sdf_mat``.

    Notes
    -----
    MATLAB stores arrays in column-major (Fortran) order; h5py reads them
    in row-major order, so a ``.T`` transpose is applied to restore the
    expected (n_trials, n_timepoints) layout.
    """
    try:
        import h5py
    except ImportError as exc:
        raise ImportError(
            "h5py is required for v7.3 .mat files: pip install h5py"
        ) from exc

    with h5py.File(str(path), 'r') as f:
        sdf_grp = f['sdf']
        return {
            'sdf': {
                'stim_on':  sdf_grp['stim_on'][:].T,
                'fixation': sdf_grp['fixation'][:].T,
                'reward':   sdf_grp['reward'][:].T,
                'brk_fix':  sdf_grp['brk_fix'][:].T,
            },
        }


def load_beh_mat(path: Path | str) -> dict:
    """
    Load a raw behavioural ``_spk.mat`` file and return trial event times.

    The file is expected to contain a MATLAB struct ``t_evt`` with fields
    corresponding to behavioural event timestamps (one value per trial,
    in seconds).

    Parameters
    ----------
    path : Path or str
        Full path to the ``_spk.mat`` file (e.g. ``spk/session_spk.mat``).

    Returns
    -------
    dict with keys:
        ``'stim_on'``      : np.ndarray, shape (n_trials,)
        ``'fixcross_fix'`` : np.ndarray, shape (n_trials,)
        ``'rew_on'``       : np.ndarray, shape (n_trials,)
        ``'brk_fix'``      : np.ndarray, shape (n_trials,)
    """
    mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    t   = mat['t_evt']

    return {
        'stim_on':      np.atleast_1d(t.stim_on),
        'fixcross_fix': np.atleast_1d(t.fixcross_fix),
        'rew_on':       np.atleast_1d(t.rew_on),
        'brk_fix':      np.atleast_1d(t.brk_fix),
    }


def load_raster_mat(path: Path | str, event: str = 'stim_on') -> list[np.ndarray]:
    """
    Load spike times from the raster field of a proc .mat file.

    Parameters
    ----------
    path : Path or str
        Full path to the proc .mat file.
    event : str
        Field name within the raster struct (e.g. 'stim_on', 'fixation').

    Returns
    -------
    list of np.ndarray
        One array per trial, each containing spike times in ms relative
        to the alignment event. Empty trials are represented as empty arrays.
    """
    mat = sio.loadmat(str(path), squeeze_me=True, struct_as_record=False)
    raster_struct = mat['raster']
    field = getattr(raster_struct, event)

    # field is a numpy object array of shape (n_trials,)
    # each element is either a 1-D array of spike times or a scalar 0 (empty)
    trials = []
    for item in np.atleast_1d(field):
        arr = np.atleast_1d(item).astype(float)
        # MATLAB empty cells come back as scalar 0 or empty array
        if arr.size == 1 and arr[0] == 0:
            trials.append(np.array([]))
        else:
            trials.append(arr)
    return trials


# ── Windows and constants ────────────────────────────────────────────────────
FF_WIN_MS      = (-600, 0)        # pre-stimulus spike count window (ms)
FIX_SDF_ZERO   = 2000             # index of t=0 in fixation-aligned SDF
FIX_BASELINE   = slice(0, 1000)   # -2000 to -1000 ms in fixation SDF
AC_LAG_MAX_MS  = 800              # maximum lag for autocorrelation fit
MIN_TRIALS     = 5                # minimum trials for reliable FF/CV
MIN_MEAN_COUNT = 0.5              # minimum mean spike count (spikes/trial)
TAU_BOUNDS     = (5, 2000)        # plausible timescale range in ms


def exp_decay(lag, A, tau_c, B):
    """Exponential decay model: A * exp(-lag / tau_c) + B."""
    return A * np.exp(-lag / tau_c) + B


def compute_ff_cv(spike_times_list: list, win: tuple, scale: float = 1.0) -> tuple:
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

def compute_timescale(sdf_fixation: np.ndarray) -> tuple[float, float, float]:
    """
    Fit intrinsic timescale from trial-averaged fixation SDF autocorrelation.

    Parameters
    ----------
    sdf_fixation : np.ndarray, shape (n_trials, n_timepoints)
        SDF aligned to fixation cross, spanning -2000 to 2000 ms.

    Returns
    -------
    tau_c : float
        Fitted timescale in ms. NaN if fit failed.
    A : float
        Fitted amplitude. NaN if fit failed.
    B : float
        Fitted offset. NaN if fit failed.
    """
    # Extract baseline window: -2000 to -1000 ms (indices 0:1000)
    baseline = np.nanmean(sdf_fixation[:, FIX_BASELINE], axis=0)  # (1000,)

    if np.all(np.isnan(baseline)) or baseline.std() < 1e-6:
        return np.nan, np.nan, np.nan

    # Normalised autocorrelation
    b      = baseline - baseline.mean()
    ac_full = np.correlate(b, b, mode='full')
    ac_full /= ac_full[len(b) - 1]          # normalise by zero-lag
    lags    = np.arange(len(b))              # 0 to 999 ms

    # Restrict to 0-800 ms
    fit_mask = lags <= AC_LAG_MAX_MS
    lags_fit = lags[fit_mask].astype(float)
    ac_fit   = ac_full[len(b) - 1:][fit_mask]   # positive lags only

    try:
        p0     = (0.5, 100.0, 0.0)
        bounds = ([0, TAU_BOUNDS[0], -1], [2, TAU_BOUNDS[1], 1])
        popt, _ = curve_fit(
            exp_decay, lags_fit, ac_fit,
            p0=p0, bounds=bounds, maxfev=5000
        )
        A_fit, tau_c_fit, B_fit = popt

        # Reject if tau_c outside plausible range (redundant with bounds
        # but kept as explicit quality gate)
        if not (TAU_BOUNDS[0] <= tau_c_fit <= TAU_BOUNDS[1]):
            return np.nan, np.nan, np.nan

        return tau_c_fit, A_fit, B_fit

    except (RuntimeError, ValueError):
        return np.nan, np.nan, np.nan


def load_mat_hdf5(path):
    """
    Load a v7.3 .mat file using h5py.
    Returns a dict with 'sdf' and 'raster' unpacked.
    """
    with h5py.File(str(path), 'r') as f:

        # ── SDF ─────────────────────────────────────────────────────────
        sdf = {
            field: f['sdf'][field][:].T.astype(float)
            for field in ('stim_on', 'fixation', 'reward', 'brk_fix')
            if field in f['sdf']
        }

        # ── Raster (cell array of object references) ─────────────────────
        # Each event field is an (n_trials, 1) array of HDF5 object refs.
        # Dereference each ref to get the spike time vector for that trial.
        raster = {}
        for event in ('stim_on', 'fixation', 'reward', 'brk_fix'):
            if event not in f['raster']:
                continue
            refs = f['raster'][event]  # shape (n_trials, 1) or (1, n_trials)
            refs = refs[:].flatten()
            trials = []
            for ref in refs:
                spk = f[ref][()].flatten().astype(float)
                trials.append(spk if spk.size > 0 else np.array([]))
            raster[event] = trials

    return {'sdf': sdf, 'raster': raster}