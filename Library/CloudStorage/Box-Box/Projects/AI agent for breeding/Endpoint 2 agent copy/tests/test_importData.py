"""
Test importData functionality
Translated from AlphaSimR-master/tests/testthat/test-importData.R
"""

import numpy as np
import pandas as pd
from AlphaSimPy import SimParam, new_pop, import_haplo


def test_import_trait():
    """Exercise import_haplo + SimParam.import_trait."""
    haplo = np.array(
        [[1, 1, 0, 1, 0], [1, 1, 0, 1, 0], [0, 1, 1, 0, 0], [0, 1, 1, 0, 0]]
    )
    gen_map = pd.DataFrame(
        {
            "markerName": ["a", "b", "c", "d", "e"],
            "chromosome": [1, 1, 1, 2, 2],
            "position": [0, 0.5, 1, 0.15, 0.4],
        }
    )
    ped = pd.DataFrame({"id": ["a", "b"], "mother": [0, 0], "father": [0, 0]})

    founder_pop = import_haplo(haplo=haplo, gen_map=gen_map, ploidy=2, ped=ped)
    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1

    my_trait = pd.DataFrame({"marker": ["a", "c", "d"], "a": [1, -1, 1]})

    SP.import_trait(
        marker_names=my_trait["marker"].values,
        add_eff=my_trait["a"].values,
        name="myTrait",
    )

    np.testing.assert_allclose(SP.traits[0].add_eff, my_trait["a"].values, atol=1e-6)
    np.testing.assert_allclose(SP.traits[0].intercept, 0, atol=1e-6)

    pop = new_pop(founder_pop, sim_param=SP)
    np.testing.assert_allclose(pop.gv[0, 0], 3, atol=1e-6)
    np.testing.assert_allclose(pop.gv[1, 0], -3, atol=1e-6)
