import numpy as np

import lalpulsar


def mean_anomaly_from_true(true_anomaly, eccentricity):
    return (
        np.arctan2(
            -np.sqrt(1 - eccentricity**2) * np.sin(true_anomaly),
            -eccentricity - np.cos(true_anomaly),
        )
        + np.pi
        - eccentricity
        * np.sqrt(1 - eccentricity**2)
        * np.sin(true_anomaly)
        / (1 + eccentricity * np.cos(true_anomaly))
    )


def eccentric_anomaly_from_mean(mean_anomaly, eccentricity):
    return lalpulsar.ComputeEccentricAnomaly(mean_anomaly, eccentricity)


def true_anomaly_from_eccentric(eccentric_anomaly, eccentricity):
    # formula from https://ui.adsabs.harvard.edu/abs/1973CeMec...7..388B/abstract
    # avoids numerical issues

    beta = eccentricity / (1 + np.sqrt(1 - eccentricity**2))
    return eccentric_anomaly + 2 * np.arctan2(
        beta * np.sin(eccentric_anomaly), 1 - beta * np.cos(eccentric_anomaly)
    )
