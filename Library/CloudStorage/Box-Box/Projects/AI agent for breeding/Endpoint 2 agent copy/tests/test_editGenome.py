"""
Test edit_genome functionality
Translated from AlphaSimR-master/tests/testthat/test-editGenome.R
"""

import numpy as np
from AlphaSimPy import edit_genome, new_map_pop, new_pop, rand_cross, run_macs, SimParam, var_g


def test_edit_genome():
    """Test edit_genome on a deterministic tiny MapPop."""
    try:
        founder_pop = new_map_pop(
            gen_map=[np.array([0])],
            haplotypes=[np.array([[0], [0], [1], [1]])],
        )
    except Exception:
        founder_pop = run_macs(n_ind=2, n_chr=1, seg_sites=1, inbred=True)

    SP = SimParam(founder_pop=founder_pop)
    SP.n_threads = 1
    SP.add_trait_a(n_qtl_per_chr=1, mean=0, var=1)
    pop = new_pop(founder_pop, sim_param=SP)
    pop = edit_genome(pop=pop, ind=0, chr=1, seg_sites=1, allele=1, sim_param=SP)
    np.testing.assert_allclose(var_g(pop), 0, atol=1e-6)
    pop = rand_cross(pop=pop, n_crosses=10, sim_param=SP)
    np.testing.assert_allclose(var_g(pop), 0, atol=1e-6)
