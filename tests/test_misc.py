"""
Test misc functionality
Translated from AlphaSimR-master/tests/testthat/test-misc.R
"""

import pytest
import numpy as np
from AlphaSimPy import merge_pops, new_multi_pop, new_pop, pull_seg_site_geno, quick_haplo, run_macs, SimParam


@pytest.mark.parametrize("quick", [False, True])
def test_misc_merge_pops(quick):
    """Misc slot handling and merge_pops stacking."""
    if quick:
        try:
            founder_pop = quick_haplo(n_ind=2, n_chr=1, seg_sites=10)
        except Exception:
            pytest.skip("quick_haplo unavailable")
    else:
        try:
            founder_pop = run_macs(n_ind=2, n_chr=1, seg_sites=10, inbred=True)
        except RuntimeError:
            founder_pop = quick_haplo(n_ind=2, n_chr=1, seg_sites=10)
    SP = SimParam(founder_pop=founder_pop)
    SP.add_trait_a(10)
    pop_orig = new_pop(founder_pop, sim_param=SP)
    multi_pop = new_multi_pop(pop_orig, pop_orig)

    assert pop_orig.misc == {}
    assert pop_orig.misc_pop == {}

    pop_sub = pop_orig[0]
    assert pop_sub.misc == {}
    assert pop_sub.misc_pop == {}

    pop = pop_orig
    pop.misc["vec"] = np.random.normal(size=2)
    pop.misc["mat"] = np.array([[1, 2], [3, 4]])
    pop.misc["counts"] = np.array([float(pop_orig.n_ind), float(len(multi_pop))])
    pop.misc_pop["vec"] = np.sum(pop.misc["vec"])
    pop.misc_pop["af"] = np.mean(pull_seg_site_geno(pop, sim_param=SP), axis=0)
    pop.misc_pop["mat"] = np.array([[1, 2], [3, 4]])

    pop_c = merge_pops([pop, pop])
    np.testing.assert_array_equal(pop_c.misc["vec"], np.concatenate([pop.misc["vec"], pop.misc["vec"]]))
    np.testing.assert_array_equal(pop_c.misc["mat"], np.vstack([pop.misc["mat"], pop.misc["mat"]]))
    assert pop_c.misc_pop == {}

    pop_c = merge_pops([pop, pop[[0]]])
    np.testing.assert_array_equal(pop_c.misc["vec"], np.concatenate([pop.misc["vec"], pop.misc["vec"][0:1]]))
    np.testing.assert_array_equal(pop_c.misc["mat"], np.vstack([pop.misc["mat"], pop.misc["mat"][0:1, :]]))
    assert pop_c.misc_pop == {}
