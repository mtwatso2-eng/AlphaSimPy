"""
Test statistics functionality
Translated from AlphaSimR-master/tests/testthat/test-statistics.R
"""

import pytest
import numpy as np
from AlphaSimPy import add_error

# These tests will fail with a small probability and/or 
# they are slower tests, so they are skipped on CRAN.
# If any of these test fail, rerun the tests multiple times.
# Frequent failures indicate a problem.


@pytest.mark.skip(reason="Skipped on CRAN - slow/stochastic test")
def test_add_error():
    """Test add_error function"""
    gv = np.zeros((10000, 2))
    var_e = np.array([1, 1])
    pheno = add_error(gv, var_e, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e), rtol=0.1)

    var_e = np.eye(2)
    pheno = add_error(gv, var_e, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e), rtol=0.1)

    pheno = add_error(gv, var_e, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e) / 4, rtol=0.1)

    var_e = 0.5 * np.eye(2) + 0.5
    pheno = add_error(gv, var_e, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e), rtol=0.1)

    pheno = add_error(gv, var_e, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e) / 4, rtol=0.1)

    var_e = 1.5 * np.eye(2) - 0.5
    pheno = add_error(gv, var_e, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e), rtol=0.1)

    pheno = add_error(gv, var_e, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(var_e) / 4, rtol=0.1)

    var_e = np.array([[1, 0.5], [-0.5, 1]])
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        add_error(gv, var_e, reps=1)

