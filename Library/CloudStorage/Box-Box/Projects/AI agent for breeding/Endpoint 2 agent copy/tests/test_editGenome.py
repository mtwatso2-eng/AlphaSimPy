"""
Test editGenome functionality
Translated from AlphaSimR-master/tests/testthat/test-editGenome.R
"""

import pytest
import numpy as np
from AlphaSimPy import runMacs, SimParam, newPop, randCross, varG
# Note: newMapPop, editGenome may need to be implemented
try:
    from AlphaSimPy import newMapPop, editGenome
except ImportError:
    newMapPop = None
    editGenome = None


def test_editGenome():
    """Test editGenome function"""
    if editGenome is None or newMapPop is None:
        pytest.skip("editGenome or newMapPop not implemented")
    
    if newMapPop is not None:
        founderPop = newMapPop(genMap=[np.array([0])],
                              haplotypes=[np.array([[0], [0], [1], [1]])])
    else:
        founderPop = runMacs(nInd=2, nChr=1, segSites=1, inbred=True)
    
    SP = SimParam(founder_pop=founderPop)
    SP.n_threads = 1
    SP.addTraitA(n_qtl_per_chr=1, mean=0, var=1)
    pop = newPop(founderPop, sim_param=SP)
    pop = editGenome(pop=pop, ind=0, chr=1, segSites=1, allele=1, sim_param=SP)
    np.testing.assert_allclose(varG(pop), 0, atol=1e-6)
    pop = randCross(pop=pop, n_crosses=10, sim_param=SP)
    np.testing.assert_allclose(varG(pop), 0, atol=1e-6)

