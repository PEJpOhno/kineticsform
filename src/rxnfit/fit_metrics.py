# Copyright (c) 2025 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# license that can be found in the LICENSE file.

"""Goodness-of-fit metrics shared by solvers and fitters.

Computes RSS, TSS, R², RMSE, and MAE. Residuals and totals omit NaN pairs so
RSS, SAE, and ``n_datapoints`` always refer to the same valid (time, species)
samples.
"""

import numpy as np

from .expdata_reader import time_course, align_expdata_to_function_names


# Threshold below which TSS is considered nearly zero; callers may warn.
TSS_MIN_THRESHOLD = 1e-12


def expdata_df_to_datasets(expdata_df, function_names):
    """Convert a single experimental DataFrame to datasets format.

    Args:
        expdata_df (pandas.DataFrame): Time in first column, species
            concentrations in remaining columns. Column names must include
            all function_names.
        function_names (list[str]): Chemical species names (e.g. ODE order).

    Returns:
        list: One-element list of dicts with keys 't_list', 'C_exp_list'
            (aligned to function_names order; valid points only, NaN removed).

    Raises:
        ValueError: If a species in function_names is not in the DataFrame.
    """
    t_list, C_exp_list = time_course(expdata_df)
    columns = list(expdata_df.columns[1:])
    t_aligned, C_aligned = align_expdata_to_function_names(
        t_list, C_exp_list, columns, function_names
    )
    return [{"t_list": t_aligned, "C_exp_list": C_aligned}]


def _compute_tss(datasets):
    """Compute TSS using per-species means over all valid observations.

    Same valid points as used for RSS (each dataset's t_list[i], C_exp_list[i]).
    """
    if not datasets:
        return 0.0
    n_species = len(datasets[0]["C_exp_list"])
    all_C = [[] for _ in range(n_species)]
    for ds in datasets:
        for i in range(n_species):
            all_C[i].append(ds["C_exp_list"][i])
    # Flatten per species
    concat = []
    for i in range(n_species):
        concat.append(np.concatenate([np.ravel(c) for c in all_C[i]]))
    tss = 0.0
    for i in range(n_species):
        arr = concat[i]
        if arr.size == 0:
            continue
        mu = np.mean(arr)
        tss += np.sum((arr - mu) ** 2)
    return float(tss)


def fit_metrics(datasets, rss, sae, n_datapoints):
    """Compute fit metrics from datasets and residual sums.

    TSS uses per-species means over the same valid points as RSS/SAE.
    Callers should pass keyword arguments for ``sae`` and ``n_datapoints``::

        fit_metrics(datasets, rss, sae=sae, n_datapoints=n)

    Args:
        datasets (list[dict]): List of dicts with 't_list', 'C_exp_list'
            (each list of arrays in species order; valid points only).
        rss (float): Residual sum of squares (already computed).
        sae (float): Sum of absolute residuals (already computed).
        n_datapoints (int): Number of valid residual points (same set as RSS/SAE).

    Returns:
        dict: Keys 'rss', 'tss', 'r2', 'rmse', 'mae', 'n_datapoints'.
            If TSS <= 0, 'r2' is NaN. If ``n_datapoints`` <= 0, 'rmse' and
            'mae' are NaN. Callers should warn when tss < TSS_MIN_THRESHOLD.
    """
    rss_f = float(rss)
    sae_f = float(sae)
    n = int(n_datapoints)
    tss = _compute_tss(datasets)
    if tss <= 0:
        r2 = np.nan
    else:
        r2 = 1.0 - (rss_f / tss)
    if n <= 0:
        rmse = np.nan
        mae = np.nan
    else:
        rmse = float(np.sqrt(rss_f / n))
        mae = float(sae_f / n)
    return {
        "rss": rss_f,
        "tss": float(tss),
        "r2": float(r2),
        "rmse": rmse,
        "mae": mae,
        "n_datapoints": n,
    }
