"""
Test phenotypes functionality
Translated from AlphaSimR-master/tests/testthat/test-phenotypes.R
"""

import pytest
import numpy as np
# Note: asCategorical may need to be implemented
try:
    from AlphaSimPy import asCategorical
except ImportError:
    # If not implemented, skip tests
    pytest.skip("asCategorical not implemented", allow_module_level=True)


def test_asCategorical_converts_correctly():
    """Test asCategorical function converts correctly"""
    cont = np.zeros((7, 3))
    cont[:, 0] = np.array([-3, -2, -1, 0, 1, 2, 3])
    cont[:, 1] = np.array([-3, -2, -1, 0, 1, 2, 3])
    cont[:, 2] = np.array([-3, -2, -1, 0, 1, 2, 3])

    result1 = asCategorical(x=cont[:, 0])
    expected1 = np.array([[1], [1], [1], [2], [2], [2], [2]])
    np.testing.assert_array_equal(result1, expected1)

    result2 = asCategorical(x=cont[:, 0], threshold=np.array([-1, 0, 1]))
    expected2 = np.array([[np.nan], [np.nan], [1], [2], [2], [np.nan], [np.nan]])
    # equal_nan added in numpy 1.19; use manual check for older versions
    assert (np.isnan(result2) == np.isnan(expected2)).all()
    np.testing.assert_array_equal(np.nan_to_num(result2), np.nan_to_num(expected2))

    result3 = asCategorical(x=cont[:, 0], threshold=np.array([-np.inf, -1, 0, 1, np.inf]))
    expected3 = np.array([[1], [1], [2], [3], [4], [4], [4]])
    np.testing.assert_array_equal(result3, expected3)

    # p=0.5 (scalar) is auto-converted to [0.5, 0.5]
    result4 = asCategorical(x=cont[:, 0], p=0.5)
    result4_expected = asCategorical(x=cont[:, 0], p=np.array([0.5, 0.5]))
    np.testing.assert_array_equal(result4, result4_expected)

    trt_mean = np.mean(cont, axis=0)
    trt_var = np.var(cont, axis=0)
    
    result5 = asCategorical(x=cont[:, 0], p=np.array([0.5, 0.5]), var=trt_var[0])
    expected5 = np.array([[1], [1], [1], [2], [2], [2], [2]])
    np.testing.assert_array_equal(result5, expected5)

    result6 = asCategorical(x=cont[:, 0], p=np.array([2/7, 1/7, 1/7, 3/7]), var=trt_var[0])
    expected6 = np.array([[1], [1], [2], [3], [4], [4], [4]])
    np.testing.assert_array_equal(result6, expected6)

    # Test error for matrix input
    with pytest.raises((ValueError, TypeError)):
        asCategorical(x=cont)

    cont2 = asCategorical(x=cont,
                          threshold=[None,
                                     np.array([-np.inf, 0, np.inf]),
                                     None])
    cont2_exp = cont.copy()
    cont2_exp[:, 1] = np.array([1, 1, 1, 2, 2, 2, 2])
    np.testing.assert_array_equal(cont2, cont2_exp)

    # Test error for p without mean/var
    with pytest.raises((ValueError, TypeError)):
        asCategorical(x=cont, p=np.array([0.5, 0.5]))
    
    p_list = [None, np.array([0.5, 0.5]), None]
    with pytest.raises((ValueError, TypeError)):
        asCategorical(x=cont, p=p_list)
    
    with pytest.raises((ValueError, TypeError)):
        asCategorical(x=cont, p=p_list, mean=trt_mean)
    
    cont2 = asCategorical(x=cont, p=p_list, mean=trt_mean, var=trt_var)
    cont2_exp = cont.copy()
    cont2_exp[:, 1] = np.array([1, 1, 1, 2, 2, 2, 2])
    np.testing.assert_array_equal(cont2, cont2_exp)

