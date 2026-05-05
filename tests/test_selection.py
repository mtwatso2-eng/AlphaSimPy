"""
Test selection functionality
Translated from AlphaSimR-master/tests/testthat/test-selection.R
"""

import pytest
import numpy as np
from AlphaSimPy import quick_haplo, run_macs, SimParam, new_pop, select_ind


def test_select_ind_and_get_response():
    """Test select_ind and getResponse functions"""
    try:
        founder_pop = quick_haplo(n_ind=10, n_chr=1, seg_sites=10)
    except Exception:
        founder_pop = run_macs(n_ind=10, n_chr=1, seg_sites=10, inbred=True)

    SP = SimParam(founder_pop=founder_pop)
    SP.add_trait_a(10)
    SP.setVarE(h2=0.5)
    pop = new_pop(founder_pop, sim_param=SP)

    pop2 = select_ind(pop, 5, sim_param=SP)
    sorted_indices = np.argsort(pop.pheno[:, 0])[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop2.id == expected_ids

    pop2b = select_ind(pop, 5, trait="Trait1", sim_param=SP)
    assert pop2.id == pop2b.id

    def squared_deviation(x, optima=0):
        return (x - optima) ** 2

    pop3 = select_ind(pop, 5, trait=squared_deviation, select_top=True, sim_param=SP)
    sorted_indices = np.argsort(squared_deviation(pop.pheno[:, 0]))[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop3.id == expected_ids

    pop4 = select_ind(pop, 5, trait=squared_deviation, select_top=False, sim_param=SP)
    sorted_indices = np.argsort(squared_deviation(pop.pheno[:, 0]))[:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop4.id == expected_ids

    pop.misc = {"smth": np.random.normal(size=10), "smth2": np.random.normal(size=10)}

    def use_func(pop, trait=None):
        return pop.misc["smth"] + pop.misc["smth2"]

    pop5 = select_ind(pop, 5, use=use_func, sim_param=SP)
    sorted_indices = np.argsort(use_func(pop))[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop5.id == expected_ids

    def use_func2(pop, trait=None):
        return np.column_stack([pop.misc["smth"], pop.misc["smth2"]])

    def trt_func(x):
        return np.sum(x, axis=1)

    pop6 = select_ind(pop, 5, trait=trt_func, use=use_func2, sim_param=SP)
    assert pop5.id == pop6.id
