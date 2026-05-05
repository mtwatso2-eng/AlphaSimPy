"""
Test hybrids functionality
Translated from AlphaSimR-master/tests/testthat/test-hybrids.R
"""

import numpy as np
from AlphaSimPy import calc_gca, hybrid_cross, new_pop, quick_haplo, run_macs, select_ind, SimParam


def _founder_two():
    try:
        return run_macs(n_ind=2, n_chr=1, seg_sites=1, inbred=True)
    except RuntimeError:
        return quick_haplo(n_ind=2, n_chr=1, seg_sites=1, inbred=True)


def test_hybrid_cross():
    """Test hybrid_cross function"""
    founder_pop = _founder_two()

    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    pop = new_pop(founder_pop, sim_param=SP)

    hybrid = hybrid_cross(pop, pop, sim_param=SP)
    assert hybrid.n_ind == 4

    hybrid = hybrid_cross(pop, pop, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 4

    pop_single = select_ind(pop, 1, sim_param=SP)
    hybrid = hybrid_cross(pop, pop_single, sim_param=SP)
    assert hybrid.n_ind == 2

    hybrid = hybrid_cross(pop, pop_single, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 2

    hybrid = hybrid_cross(pop_single, pop_single, sim_param=SP)
    assert hybrid.n_ind == 1

    hybrid = hybrid_cross(pop_single, pop_single, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 1


def test_calc_gca():
    """Test calc_gca function"""
    founder_pop = _founder_two()

    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=np.array([0, 0]), var=np.array([1, 1]))
    SP.setVarE(varE=np.array([1, 1]))
    pop = new_pop(founder_pop, sim_param=SP)

    hybrid = hybrid_cross(pop, pop, return_hybrid_pop=True, sim_param=SP)
    gca = calc_gca(hybrid)
    assert gca["GCAf"].shape[0] == 2
    assert gca["GCAm"].shape[0] == 2
    assert gca["SCA"].shape[0] == 4

    pop_single = select_ind(pop, 1, sim_param=SP)
    hybrid = hybrid_cross(pop, pop_single, return_hybrid_pop=True, sim_param=SP)
    gca = calc_gca(hybrid)
    assert gca["GCAf"].shape[0] == 2
    assert gca["GCAm"].shape[0] == 1
    assert gca["SCA"].shape[0] == 2

    hybrid = hybrid_cross(pop_single, pop, return_hybrid_pop=True, sim_param=SP)
    gca = calc_gca(hybrid)
    assert gca["GCAf"].shape[0] == 1
    assert gca["GCAm"].shape[0] == 2
    assert gca["SCA"].shape[0] == 2

    hybrid = hybrid_cross(pop_single, pop_single, return_hybrid_pop=True, sim_param=SP)
    gca = calc_gca(hybrid)
    assert gca["GCAf"].shape[0] == 1
    assert gca["GCAm"].shape[0] == 1
    assert gca["SCA"].shape[0] == 1
