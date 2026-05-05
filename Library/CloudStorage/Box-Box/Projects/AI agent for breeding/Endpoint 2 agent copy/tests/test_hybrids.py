"""
Test hybrids functionality
Translated from AlphaSimR-master/tests/testthat/test-hybrids.R
"""

import pytest
import numpy as np
# Note: newMapPop, quickHaplo may need to be implemented
# Using runMacs as alternative for now
from AlphaSimPy import runMacs, SimParam, newPop, hybridCross, calcGCA, selectInd


def test_hybridCross():
    """Test hybridCross function"""
    # Create founder population with 2 individuals, 1 chromosome, 1 seg site
    # Using runMacs as newMapPop may not be implemented
    founderPop = runMacs(nInd=2, nChr=1, segSites=1, inbred=True)
    
    SP = SimParam(founder_pop=founderPop)
    SP.n_threads = 1
    SP.addTraitA(n_qtl_per_chr=1, mean=0, var=1)
    pop = newPop(founderPop, sim_param=SP)
    
    # 2x2
    hybrid = hybridCross(pop, pop, sim_param=SP)
    assert hybrid.n_ind == 4
    
    hybrid = hybridCross(pop, pop, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 4
    
    # 2x1 - Create single-individual population using selectInd
    pop_single = selectInd(pop, 1, sim_param=SP)
    hybrid = hybridCross(pop, pop_single, sim_param=SP)
    assert hybrid.n_ind == 2
    
    hybrid = hybridCross(pop, pop_single, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 2
    
    # 1x1
    hybrid = hybridCross(pop_single, pop_single, sim_param=SP)
    assert hybrid.n_ind == 1
    
    hybrid = hybridCross(pop_single, pop_single, return_hybrid_pop=True, sim_param=SP)
    assert hybrid.n_ind == 1


def test_calcGCA():
    """Test calcGCA function"""
    founderPop = runMacs(nInd=2, nChr=1, segSites=1, inbred=True)
    
    SP = SimParam(founder_pop=founderPop)
    SP.n_threads = 1
    SP.addTraitA(n_qtl_per_chr=1, mean=np.array([0, 0]), var=np.array([1, 1]))
    SP.setVarE(var_e=np.array([1, 1]))
    pop = newPop(founderPop, sim_param=SP)
    
    # 2x2
    hybrid = hybridCross(pop, pop, return_hybrid_pop=True, sim_param=SP)
    GCA = calcGCA(hybrid)
    assert GCA['GCAf'].shape[0] == 2
    assert GCA['GCAm'].shape[0] == 2
    assert GCA['SCA'].shape[0] == 4
    
    # 2x1 - Create single-individual population using selectInd
    pop_single = selectInd(pop, 1, sim_param=SP)
    hybrid = hybridCross(pop, pop_single, return_hybrid_pop=True, sim_param=SP)
    GCA = calcGCA(hybrid)
    assert GCA['GCAf'].shape[0] == 2
    assert GCA['GCAm'].shape[0] == 1
    assert GCA['SCA'].shape[0] == 2
    
    # 1x2
    hybrid = hybridCross(pop_single, pop, return_hybrid_pop=True, sim_param=SP)
    GCA = calcGCA(hybrid)
    assert GCA['GCAf'].shape[0] == 1
    assert GCA['GCAm'].shape[0] == 2
    assert GCA['SCA'].shape[0] == 2
    
    # 1x1
    hybrid = hybridCross(pop_single, pop_single, return_hybrid_pop=True, sim_param=SP)
    GCA = calcGCA(hybrid)
    assert GCA['GCAf'].shape[0] == 1
    assert GCA['GCAm'].shape[0] == 1
    assert GCA['SCA'].shape[0] == 1

