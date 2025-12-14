import pytest
from ..models import teobresums
from ..models import pyseobnr_model
import astropy.units as u

ALL_MODELS = {
    "TEOBResumSDALI": teobresums.TEOBResumSDALI,
    "SEOBNRv5HM": pyseobnr_model.SEOBNRv5HM,
}


def pytest_addoption(parser):
    parser.addoption(
        "--plot",
        action="store_true",
        default=False,
        help="""Produce plots for the tests that support them.
        """,
    )

    available_models_str = ", ".join(list(ALL_MODELS.keys()))

    parser.addoption(
        "--model",
        action="append",
        default=[],
        help=f"""
        Model to test.
        If no model is specified, all available models are tested.
        Options: {available_models_str}.
        It is possible to specify multiple models, e.g. --model=TEOBResumSDALI --model=SEOBNRv5HM""",
    )


@pytest.fixture
def plot(request):
    return request.config.getoption("--plot")


def pytest_generate_tests(metafunc):
    if "gen" in metafunc.fixturenames:
        models = metafunc.config.getoption("model")

        if len(models) == 0:
            models_to_test = ALL_MODELS
        else:
            try:
                models_to_test = {model: ALL_MODELS[model] for model in models}
            except KeyError as e:
                raise ValueError("Model not found") from e
        metafunc.parametrize(
            "gen", models_to_test.values(), indirect=True, ids=models_to_test.keys()
        )


@pytest.fixture(scope="module")
def gen(request):
    # this is the generator to be used in all tests;
    # it is determined by the pytest_generate_tests function
    return request.param()


@pytest.fixture
def parameters():
    return {
        "mass1": 60.0 * u.solMass,
        "mass2": 50.0 * u.solMass,
        "spin1x": 0.0 * u.dimensionless_unscaled,
        "spin1y": 0.0 * u.dimensionless_unscaled,
        "spin1z": -0.1 * u.dimensionless_unscaled,
        "spin2x": 0.0 * u.dimensionless_unscaled,
        "spin2y": 0.0 * u.dimensionless_unscaled,
        "spin2z": 0.1 * u.dimensionless_unscaled,
        "eccentricity": 0.5 * u.dimensionless_unscaled,
        "deltaT": 1.0 / 1024.0 / 16 * u.s,
        "f22_start": 20.0 * u.Hz,
        "f22_ref": 20.0 * u.Hz,
        "distance": 1000.0 * u.Mpc,
        "inclination": 1.0 * u.rad,
        "phi_ref": -1.0 * u.rad,
        "longAscNodes": 0.0 * u.rad,
        "meanPerAno": 0.0 * u.rad,
        "condition": 0,
    }
