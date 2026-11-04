"""Display switch of the Python tests, as TEST_PLOT_ON / PLOT_ON for the C++ tests.

The tests draw nothing by default: an interactive window blocks the suite until it is closed by
hand. Set the environment variable GBS_TEST_PLOT_ON (to anything but "", "0" or "false") to see
the plots.
"""
import os

PLOT_ON = os.environ.get("GBS_TEST_PLOT_ON", "").strip().lower() not in ("", "0", "false")
