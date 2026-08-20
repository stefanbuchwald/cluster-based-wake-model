import unittest

import numpy as np

from rom.pod import build_component_weights, weighted_snapshot_pod


class WeightedPodTests(unittest.TestCase):
    def test_component_weights_repeat_cell_volumes_per_component(self):
        volumes = np.array([1.0, 4.0])
        layout = {
            "fields": [
                {"name": "U", "n_components": 3},
                {"name": "p", "n_components": 1},
            ],
            "n_features": 8,
        }
        weights = build_component_weights(volumes, layout)
        np.testing.assert_allclose(weights, [1, 1, 1, 4, 4, 4, 1, 4])

    def test_weighted_energy_matches_discrete_volume_integral(self):
        X = np.array([[0.0, 0.0], [2.0, 1.0], [-2.0, -1.0]])
        weights = np.array([1.0, 4.0])
        pod = weighted_snapshot_pod(X, weights=weights, energy_threshold=1.0)
        centered = X - X.mean(axis=0)
        expected_total = np.sum((centered**2) * weights[None, :]) / (X.shape[0] - 1)
        self.assertAlmostEqual(float(np.sum(pod["eigenvalues"])), expected_total)

    def test_weighted_and_unweighted_pod_differ_for_nonuniform_volumes(self):
        X = np.array([[2.0, 0.0], [0.0, 1.0], [-2.0, 0.0], [0.0, -1.0]])
        unweighted = weighted_snapshot_pod(X, weights=np.ones(2), n_modes_max=1)
        weighted = weighted_snapshot_pod(X, weights=np.array([1.0, 10.0]), n_modes_max=1)
        self.assertFalse(
            np.allclose(
                np.abs(unweighted["spatial_modes"][:, 0]),
                np.abs(weighted["spatial_modes"][:, 0]),
            )
        )


if __name__ == "__main__":
    unittest.main()
