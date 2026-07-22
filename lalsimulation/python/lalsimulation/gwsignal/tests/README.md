# Automatic pre-review tests

This suite of tests is meant to simplify the task of performing certain
pre-review checks.
From the [pre-review checklist](https://git.ligo.org/waveforms/1-main/-/wikis/pre-review-checklist), the [standard waveform checks](https://git.ligo.org/waveforms/reviews/seobnrv5/-/wikis/standard-waveform-tests) can be automated in the manner described here.

The tests are launched with `pytest`; when launching them, additional options can also be given, as shown by `pytest --help`.
A few custom options are specifically defined for this test suite, in `conftest.py`:

```
Custom options:
  --plot                Produce plots for the tests that support them.
  --model=MODEL         Model to test.
                        If no model is specified, all available models are tested.
                        Options: TEOBResumSDALI, SEOBNRv5HM.
                        It is possible to specify multiple models, e.g. --model=TEOBResumSDALI --model=SEOBNRv5HM
```

The model parametrization is achieved thanks to the `gen` fixture,
which is defined as the tests are launched --- almost every test will
use this fixture, which is defined to cycle across all models meant to be tested.
