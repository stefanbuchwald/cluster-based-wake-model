import unittest

import numpy as np

from rom.transitions import (
    cnm_matrices,
    deterministic_cnm,
    enumerate_cycles,
    probabilistic_cmm,
    probabilistic_cnm,
)


class TransitionModelTests(unittest.TestCase):
    def test_cmm_samples_p_and_advances_exact_dt(self):
        P = np.array([[0.75, 0.2], [0.25, 0.8]])
        centers = np.array([[0.0], [1.0]])
        result = probabilistic_cmm(
            P=P,
            centers=centers,
            start_cluster=0,
            n_steps=40000,
            dt=0.25,
            rng=np.random.default_rng(123),
        )
        self.assertTrue(np.allclose(np.diff(result["time"]), 0.25))
        labels = result["labels"]
        counts = np.zeros_like(P)
        for old, new in zip(labels[:-1], labels[1:]):
            counts[new, old] += 1
        empirical = counts / counts.sum(axis=0, keepdims=True)
        np.testing.assert_allclose(empirical, P, atol=0.02)

    def test_cnm_sampling_uses_q_and_transition_times(self):
        Q = np.array([[0.0, 1.0], [1.0, 0.0]])
        T = np.array([[np.nan, 2.0], [3.0, np.nan]])
        centers = np.array([[0.0], [1.0]])
        det = deterministic_cnm(Q=Q, T=T, centers=centers, start_cluster=0, n_transitions=3)
        self.assertEqual(det["labels"].tolist(), [0, 1, 0, 1])
        np.testing.assert_allclose(np.diff(det["event_time"]), [3.0, 2.0, 3.0])
        prob = probabilistic_cnm(
            Q=Q,
            T=T,
            centers=centers,
            start_cluster=0,
            n_transitions=3,
            rng=np.random.default_rng(1),
        )
        self.assertEqual(prob["labels"].tolist(), [0, 1, 0, 1])
        np.testing.assert_allclose(np.diff(prob["event_time"]), [3.0, 2.0, 3.0])

    def test_censored_endpoint_runs_are_excluded_from_transition_times(self):
        labels = np.array([0, 0, 1, 1, 1, 0, 0])
        time = np.arange(len(labels), dtype=float)
        result = cnm_matrices(labels, time, K=2)
        self.assertTrue(result["left_censored"])
        self.assertTrue(result["right_censored"])
        self.assertTrue(np.isnan(result["T"][1, 0]))
        self.assertTrue(np.isnan(result["T"][0, 1]))
        self.assertEqual(result["n_complete_runs"], 1)

    def test_cycle_frequency_from_network_times_and_support(self):
        Q = np.array([[0.0, 1.0], [1.0, 0.0]])
        T = np.array([[np.nan, 2.0], [3.0, np.nan]])
        counts = np.array([[0, 5], [7, 0]])
        cycles = enumerate_cycles(Q=Q, T=T, counts=counts, min_count=1)
        self.assertEqual(len(cycles), 1)
        cycle = cycles[0]
        self.assertEqual(cycle["nodes"], [0, 1])
        self.assertAlmostEqual(cycle["period"], 5.0)
        self.assertAlmostEqual(cycle["frequency"], 0.2)
        self.assertEqual(cycle["min_support"], 5)


if __name__ == "__main__":
    unittest.main()
