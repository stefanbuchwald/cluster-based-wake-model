import numpy as np
import pytest

from run_cnm_analysis import prepare_time_series, probabilistic_cmm_prediction


def test_uniform_series_is_preserved():
    time = 140.0 + 0.01 * np.arange(100)
    X = np.column_stack([time, 2.0 * time])
    actual_time, actual_X = prepare_time_series(time, X)
    np.testing.assert_array_equal(actual_time, time)
    np.testing.assert_array_equal(actual_X, X)


def test_nonuniform_series_requires_explicit_resampling():
    with pytest.raises(ValueError, match="--resample-dt"):
        prepare_time_series(np.array([0.0, 0.1, 0.3]), np.ones((3, 1)))


def test_resampling_interpolates_all_features_without_extending_window():
    time = np.array([2.0, 2.1, 2.35, 2.55])
    X = np.column_stack([2.0 * time, -time + 3.0])
    regular_time, regular_X = prepare_time_series(time, X, resample_dt=0.1)
    np.testing.assert_allclose(regular_time, np.arange(2.0, 2.51, 0.1))
    np.testing.assert_allclose(regular_X, np.column_stack([2.0 * regular_time, -regular_time + 3.0]))
    assert regular_time[-1] <= time[-1]


@pytest.mark.parametrize("resample_dt", [0.0, -0.1, np.nan, np.inf, 2.0])
def test_invalid_resampling_interval_is_rejected(resample_dt):
    with pytest.raises(ValueError, match="--resample-dt"):
        prepare_time_series(np.array([0.0, 0.5, 1.0]), np.ones((3, 1)), resample_dt)


def test_cmm_rejects_nonuniform_times_when_called_directly():
    with pytest.raises(ValueError, match="uniformly sampled"):
        probabilistic_cmm_prediction(
            time=np.array([0.0, 0.1, 0.3]),
            centers=np.array([[0.0], [1.0]]),
            P=np.eye(2),
            T=np.zeros((2, 2)),
            start_cluster=0,
            n_runs=1,
            rng=np.random.default_rng(1),
        )
