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
