"""
Test selection functionality
Translated from AlphaSimR-master/tests/testthat/test-selection.R
"""

import pytest
import numpy as np
from AlphaSimPy import runMacs, SimParam, newPop, selectInd
# Note: quickHaplo may need to be implemented
try:
    from AlphaSimPy import quickHaplo
except ImportError:
    quickHaplo = None


def test_selectInd_and_getResponse():
    """Test selectInd and getResponse functions"""
    if quickHaplo is not None:
        founderPop = quickHaplo(nInd=10, nChr=1, segSites=10)
    else:
        founderPop = runMacs(nInd=10, nChr=1, segSites=10, inbred=True)
    
    SP = SimParam(founderPop)
    SP.addTraitA(10)
    SP.setVarE(h2=0.5)
    pop = newPop(founderPop, sim_param=SP)

    pop2 = selectInd(pop, 5, sim_param=SP)
    sorted_indices = np.argsort(pop.pheno[:, 0])[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop2.id == expected_ids

    pop2b = selectInd(pop, 5, trait="Trait1", sim_param=SP)
    assert pop2.id == pop2b.id

    def squaredDeviation(x, optima=0):
        return (x - optima) ** 2
    
    pop3 = selectInd(pop, 5, trait=squaredDeviation, select_top=True, sim_param=SP)
    sorted_indices = np.argsort(squaredDeviation(pop.pheno[:, 0]))[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop3.id == expected_ids

    pop4 = selectInd(pop, 5, trait=squaredDeviation, select_top=False, sim_param=SP)
    sorted_indices = np.argsort(squaredDeviation(pop.pheno[:, 0]))[:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop4.id == expected_ids

    pop.misc = {'smth': np.random.normal(size=10), 'smth2': np.random.normal(size=10)}
    
    def useFunc(pop, trait=None):
        return pop.misc['smth'] + pop.misc['smth2']
    
    pop5 = selectInd(pop, 5, use=useFunc, sim_param=SP)
    sorted_indices = np.argsort(useFunc(pop))[::-1][:5]
    expected_ids = [pop.id[i] for i in sorted_indices]
    assert pop5.id == expected_ids

    def useFunc2(pop, trait=None):
        return np.column_stack([pop.misc['smth'], pop.misc['smth2']])
    
    def trtFunc(x):
        return np.sum(x, axis=1)
    
    pop6 = selectInd(pop, 5, trait=trtFunc, use=useFunc2, sim_param=SP)
    assert pop5.id == pop6.id

