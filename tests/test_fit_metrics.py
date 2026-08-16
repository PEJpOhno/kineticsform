# Copyright (c) 2025 Mitsuru Ohno
# Use of this source code is governed by a BSD-3-style
# license that can be found in the LICENSE file.

"""Tests for fit metrics and valid-point counting consistency."""

import numpy as np
import pandas as pd
import pytest

from rxnfit.fit_metrics import fit_metrics, expdata_df_to_datasets


def _n_datapoints_from_datasets(datasets, n_species):
    """Count finite (t, C) pairs the same way as expdata_fit error accumulation."""
    n = 0
    for ds in datasets:
        for i in range(n_species):
            t_i = np.asarray(ds["t_list"][i], dtype=float)
            C_i = np.asarray(ds["C_exp_list"][i], dtype=float)
            if t_i.size == 0:
                continue
            mask = np.isfinite(t_i) & np.isfinite(C_i)
            n += int(np.count_nonzero(mask))
    return n


def _n_datapoints_from_dataframe(expdata_df, function_names):
    """Count valid cells the same way as RxnODEsolver._compute_errors."""
    name_to_idx = {name: i for i, name in enumerate(function_names)}
    n = 0
    for col in expdata_df.columns[1:]:
        if col not in name_to_idx:
            continue
        c_exp = expdata_df[col].to_numpy()
        for j in range(len(expdata_df)):
            t_j = expdata_df.iloc[j, 0]
            if pd.isna(t_j) or np.isnan(c_exp[j]):
                continue
            n += 1
    return n


def test_fit_metrics_rmse_mae_and_keys():
    datasets = [
        {
            "t_list": [np.array([0.0, 1.0]), np.array([0.0, 1.0])],
            "C_exp_list": [np.array([1.0, 3.0]), np.array([2.0, 4.0])],
        }
    ]
    # Four points: means per species are 2 and 3; TSS = 2+2 = 4
    # residuals imagined: rss=2, sae=2, n=4 -> rmse=sqrt(0.5), mae=0.5, r2=1-2/4=0.5
    metrics = fit_metrics(datasets, 2.0, sae=2.0, n_datapoints=4)
    assert set(metrics) == {"rss", "tss", "r2", "rmse", "mae", "n_datapoints"}
    assert metrics["rss"] == pytest.approx(2.0)
    assert metrics["tss"] == pytest.approx(4.0)
    assert metrics["r2"] == pytest.approx(0.5)
    assert metrics["rmse"] == pytest.approx(np.sqrt(0.5))
    assert metrics["mae"] == pytest.approx(0.5)
    assert metrics["n_datapoints"] == 4


def test_fit_metrics_n_datapoints_zero_gives_nan():
    datasets = [
        {
            "t_list": [np.array([]), np.array([])],
            "C_exp_list": [np.array([]), np.array([])],
        }
    ]
    metrics = fit_metrics(datasets, 0.0, sae=0.0, n_datapoints=0)
    assert np.isnan(metrics["rmse"])
    assert np.isnan(metrics["mae"])
    assert metrics["n_datapoints"] == 0


def test_n_datapoints_consistent_between_dataset_and_dataframe_paths():
    """Same experimental valid-point set for fit (datasets) and solver (DataFrame)."""
    df = pd.DataFrame(
        {
            "time": [0.0, 1.0, 2.0, 3.0],
            "A": [1.0, np.nan, 3.0, 4.0],
            "B": [10.0, 20.0, np.nan, 40.0],
            "extra": [0.0, 0.0, 0.0, 0.0],
        }
    )
    function_names = ["A", "B"]
    datasets = expdata_df_to_datasets(df, function_names)
    n_fit = _n_datapoints_from_datasets(datasets, len(function_names))
    n_solver = _n_datapoints_from_dataframe(df, function_names)
    assert n_fit == n_solver
    # A: 3 valid, B: 3 valid
    assert n_fit == 6


def test_n_datapoints_with_nan_time_row():
    df = pd.DataFrame(
        {
            "time": [0.0, np.nan, 2.0],
            "A": [1.0, 2.0, 3.0],
            "B": [4.0, 5.0, 6.0],
        }
    )
    function_names = ["A", "B"]
    datasets = expdata_df_to_datasets(df, function_names)
    n_fit = _n_datapoints_from_datasets(datasets, len(function_names))
    n_solver = _n_datapoints_from_dataframe(df, function_names)
    assert n_fit == n_solver
