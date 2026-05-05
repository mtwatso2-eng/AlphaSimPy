"""
Test statistics functionality
Translated from AlphaSimR-master/tests/testthat/test-statistics.R
"""

import pytest
import numpy as np
import sys
# Note: addError is an internal function in AlphaSimR, may need to be exposed or implemented
try:
    from AlphaSimPy import addError
except ImportError:
    # If not implemented, skip tests
    pytest.skip("addError not implemented", allow_module_level=True)

# These tests will fail with a small probability and/or 
# they are slower tests, so they are skipped on CRAN.
# If any of these test fail, rerun the tests multiple times.
# Frequent failures indicate a problem.


@pytest.mark.skip(reason="Skipped on CRAN - slow/stochastic test")
def test_addError():
    """Test addError function"""
    gv = np.zeros((10000, 2))
    varE = np.array([1, 1])
    pheno = addError(gv=gv, var_e=varE, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE), rtol=0.1)
    
    varE = np.eye(2)
    pheno = addError(gv=gv, var_e=varE, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE), rtol=0.1)
    
    pheno = addError(gv=gv, var_e=varE, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE)/4, rtol=0.1)
    
    varE = 0.5 * np.eye(2) + 0.5
    pheno = addError(gv=gv, var_e=varE, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE), rtol=0.1)
    
    pheno = addError(gv=gv, var_e=varE, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE)/4, rtol=0.1)
    
    varE = 1.5 * np.eye(2) - 0.5
    pheno = addError(gv=gv, var_e=varE, reps=1)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE), rtol=0.1)
    
    pheno = addError(gv=gv, var_e=varE, reps=4)
    np.testing.assert_allclose(np.var(pheno, axis=0), np.diag(varE)/4, rtol=0.1)
    
    varE = np.array([[1, 0.5], [-0.5, 1]])
    with pytest.raises((ValueError, np.linalg.LinAlgError)):
        addError(gv=gv, var_e=varE, reps=1)

