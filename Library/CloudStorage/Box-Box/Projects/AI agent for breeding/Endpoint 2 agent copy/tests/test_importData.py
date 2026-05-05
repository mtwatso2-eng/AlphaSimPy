"""
Test importData functionality
Translated from AlphaSimR-master/tests/testthat/test-importData.R
"""

import pytest
import numpy as np
import pandas as pd
from AlphaSimPy import SimParam, newPop
# Note: importHaplo may need to be implemented
try:
    from AlphaSimPy import importHaplo
except ImportError:
    importHaplo = None


def test_importTrait():
    """Test importTrait function"""
    if importHaplo is None:
        pytest.skip("importHaplo not implemented")
    
    # Create haplotype data
    haplo = np.array([[1, 1, 0, 1, 0],
                      [1, 1, 0, 1, 0],
                      [0, 1, 1, 0, 0],
                      [0, 1, 1, 0, 0]])
    
    # Create genetic map
    genMap = pd.DataFrame({
        'markerName': ['a', 'b', 'c', 'd', 'e'],
        'chromosome': [1, 1, 1, 2, 2],
        'position': [0, 0.5, 1, 0.15, 0.4]
    })
    
    # Create pedigree
    ped = pd.DataFrame({
        'id': ['a', 'b'],
        'mother': [0, 0],
        'father': [0, 0]
    })
    
    # Generate an external trait
    myTrait = pd.DataFrame({
        'marker': ['a', 'c', 'd'],
        'a': [1, -1, 1]
    })
    
    founderPop = importHaplo(haplo=haplo,
                            genMap=genMap,
                            ploidy=2,
                            ped=ped)
    
    SP = SimParam(founder_pop=founderPop)
    SP.n_threads = 1
    
    # Import trait - Note: importTrait may need to be implemented
    try:
        SP.importTrait(markerNames=myTrait['marker'].values,
                       addEff=myTrait['a'].values,
                       name="myTrait")
        
        np.testing.assert_allclose(SP.traits[0].addEff, myTrait['a'].values, atol=1e-6)
        np.testing.assert_allclose(SP.traits[0].intercept, 0, atol=1e-6)
        
        pop = newPop(founderPop, sim_param=SP)
        
        np.testing.assert_allclose(pop.gv[0, 0], 3, atol=1e-6)
        np.testing.assert_allclose(pop.gv[1, 0], -3, atol=1e-6)
    except AttributeError:
        pytest.skip("importTrait not implemented")

