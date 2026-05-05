"""
Test crossing functionality
Translated from AlphaSimR-master/tests/testthat/test-crossing.R
"""

import pytest
import numpy as np
from AlphaSimPy import (
    SimParam,
    mean_g,
    make_cross,
    make_cross2,
    make_dh,
    new_map_pop,
    new_pop,
    quick_haplo,
    rand_cross,
    rand_cross2,
    run_macs,
    select_cross,
    select_op,
    self_,
)


def _founder_map_pop():
    return new_map_pop(
        gen_map=[np.array([0])],
        haplotypes=[np.array([[1], [1], [0], [0]])],
    )


def _founder_or_macs():
    try:
        return _founder_map_pop()
    except Exception:
        try:
            return run_macs(n_ind=2, n_chr=1, seg_sites=1, inbred=True)
        except RuntimeError:
            return quick_haplo(n_ind=2, n_chr=1, seg_sites=1, inbred=True)


def test_make_cross():
    """Test make_cross function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    pop = new_pop(founder_pop, sim_param=SP)
    cross_plan = np.column_stack([np.repeat(0, 10), np.repeat(1, 10)])

    pop1 = make_cross(pop=pop, cross_plan=cross_plan, sim_param=SP)
    assert SP._last_id == 12
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)

    cross_plan = np.array([[str(x), str(y)] for x, y in zip(np.repeat(0, 10), np.repeat(1, 10))])
    pop2 = make_cross(pop=pop, cross_plan=cross_plan, sim_param=SP)
    assert SP._last_id == 22
    np.testing.assert_allclose(mean_g(pop2), 0, atol=1e-6)


def test_make_cross2():
    """Test make_cross2 function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    pop = new_pop(founder_pop, sim_param=SP)
    cross_plan = np.column_stack([np.repeat(0, 10), np.repeat(1, 10)])

    pop1 = make_cross2(females=pop, males=pop, cross_plan=cross_plan, sim_param=SP)
    assert SP._last_id == 12
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)

    cross_plan = np.array([[str(x), str(y)] for x, y in zip(np.repeat(0, 10), np.repeat(1, 10))])
    pop2 = make_cross2(females=pop, males=pop, cross_plan=cross_plan, sim_param=SP)
    assert SP._last_id == 22
    np.testing.assert_allclose(mean_g(pop2), 0, atol=1e-6)


def test_rand_cross():
    """Test rand_cross function"""
    np.random.seed(4)
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    SP.set_sexes("yes_sys")
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = rand_cross(pop=pop, n_crosses=1, n_progeny=10, sim_param=SP)
    assert SP._last_id == 12
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)


def test_rand_cross2():
    """Test rand_cross2 function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    SP.set_sexes("yes_sys")
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = rand_cross2(females=pop, males=pop, n_crosses=1, n_progeny=10, sim_param=SP)
    assert SP._last_id == 12
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)


def test_self_():
    """Test self_ (self-cross) function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = self_(pop=pop, n_progeny=1, sim_param=SP)
    assert SP._last_id == 4
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)


def test_make_dh():
    """Test make_dh function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = make_dh(pop=pop, n_dh=1, sim_param=SP)
    assert SP._last_id == 4
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)


def test_select_cross():
    """Test select_cross function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    SP.set_sexes("yes_sys")
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = select_cross(pop=pop, n_female=1, n_male=1, use="rand", n_crosses=2, sim_param=SP)
    assert SP._last_id == 4
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)


def test_select_op():
    """Test select_op function"""
    founder_pop = _founder_or_macs()
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    SP.set_track_ped(True)
    pop = new_pop(founder_pop, sim_param=SP)
    pop1 = select_op(pop=pop, n_ind=1, n_seeds=2, use="gv", sim_param=SP)
    assert SP._last_id == 4
    np.testing.assert_allclose(mean_g(pop1), 0, atol=1e-6)
