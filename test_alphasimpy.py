#!/usr/bin/env python3
"""
Test script for AlphaSimPy runMacs function
"""

import sys
import numpy as np
from AlphaSimPy import runMacs, MapPop, get_num_threads

def test_basic_functionality():
    """Test basic runMacs functionality"""
    print("Testing basic runMacs functionality...")
    
    try:
        # Test 1: Basic simulation
        founder_pop = runMacs(nInd=5, nChr=1, segSites=50)
        
        print(f"✓ Basic simulation successful")
        print(f"  - Number of individuals: {founder_pop.n_ind}")
        print(f"  - Number of chromosomes: {founder_pop.n_chr}")
        print(f"  - Ploidy: {founder_pop.ploidy}")
        print(f"  - Number of loci: {founder_pop.n_loci}")
        print(f"  - Is inbred: {founder_pop.inbred}")
        print(f"  - Genetic map length: {len(founder_pop.gen_map[0])}")
        print(f"  - Genotype array shape: {founder_pop.geno[0].shape}")
        
        # Test 2: Multiple chromosomes
        founder_pop2 = runMacs(nInd=3, nChr=2, segSites=[30, 40])
        print(f"✓ Multi-chromosome simulation successful")
        print(f"  - Chromosomes: {founder_pop2.n_chr}")
        print(f"  - Loci per chromosome: {founder_pop2.n_loci}")
        
        # Test 3: Different species
        founder_pop3 = runMacs(nInd=4, nChr=1, segSites=25, species="CATTLE")
        print(f"✓ Cattle species simulation successful")
        
        # Test 4: Inbred population
        founder_pop4 = runMacs(nInd=3, nChr=1, segSites=20, inbred=True)
        print(f"✓ Inbred simulation successful")
        print(f"  - Is inbred: {founder_pop4.inbred}")
        
        return True
        
    except Exception as e:
        print(f"✗ Test failed: {e}")
        return False

def test_parameter_validation():
    """Test parameter validation"""
    print("\nTesting parameter validation...")
    
    try:
        # Test invalid species
        try:
            runMacs(nInd=2, nChr=1, segSites=10, species="INVALID")
            print("✗ Should have failed with invalid species")
            return False
        except ValueError as e:
            print(f"✓ Correctly caught invalid species: {e}")
        
        # Test manual command without gen_len
        try:
            runMacs(nInd=2, nChr=1, segSites=10, manualCommand="1E8 -t 1E-5")
            print("✗ Should have failed without manualGenLen")
            return False
        except ValueError as e:
            print(f"✓ Correctly caught missing manualGenLen: {e}")
        
        return True
        
    except Exception as e:
        print(f"✗ Parameter validation test failed: {e}")
        return False

def test_threading():
    """Test threading functionality"""
    print("\nTesting threading...")
    
    try:
        n_threads = get_num_threads()
        print(f"✓ Available threads: {n_threads}")
        
        # Test with explicit thread count
        founder_pop = runMacs(nInd=3, nChr=1, segSites=15, nThreads=1)
        print(f"✓ Single-threaded simulation successful")
        
        return True
        
    except Exception as e:
        print(f"✗ Threading test failed: {e}")
        return False

def test_map_pop_validation():
    """Test MapPop object validation"""
    print("\nTesting MapPop validation...")
    
    try:
        # This should work
        founder_pop = runMacs(nInd=2, nChr=1, segSites=10)
        print(f"✓ Valid MapPop created successfully")
        
        # Test MapPop properties
        assert founder_pop.n_ind == 2
        assert founder_pop.n_chr == 1
        assert founder_pop.ploidy == 2
        assert len(founder_pop.n_loci) == 1
        assert len(founder_pop.geno) == 1
        assert len(founder_pop.gen_map) == 1
        assert len(founder_pop.centromere) == 1
        print(f"✓ MapPop properties are correct")
        
        return True
        
    except Exception as e:
        print(f"✗ MapPop validation test failed: {e}")
        return False

def test_newpop_functionality():
    """Test newPop functionality"""
    try:
        print("Testing newPop functionality...")
        
        # Create founder population
        founder_pop = runMacs(nInd=10, nChr=1, segSites=30)
        print(f"✓ Created founder population: {founder_pop.n_ind} individuals")
        
        # Create SimParam
        from AlphaSimPy import SimParam
        SP = SimParam(founder_pop)
        SP.addTraitA(3)  # 3 QTL per chromosome
        SP.setSexes("yes_sys")
        print(f"✓ Created SimParam: {SP.n_traits} traits")
        
        # Test newPop
        from AlphaSimPy import newPop
        pop = newPop(founder_pop, sim_param=SP)
        print(f"✓ Created population: {pop.n_ind} individuals, {pop.n_traits} traits")
        
        # Test Pop properties
        assert pop.n_ind == 10
        assert pop.n_traits == 1
        assert len(pop.id) == 10
        assert len(pop.sex) == 10
        assert pop.gv.shape == (10, 1)
        assert pop.pheno.shape == (10, 1)
        print(f"✓ Pop properties are correct")
        
        return True
        
    except Exception as e:
        print(f"✗ newPop test failed: {e}")
        return False

def main():
    """Run all tests"""
    print("AlphaSimPy Test Suite")
    print("=" * 50)
    
    tests = [
        test_basic_functionality,
        test_parameter_validation,
        test_threading,
        test_map_pop_validation,
        test_newpop_functionality
    ]
    
    passed = 0
    total = len(tests)
    
    for test in tests:
        if test():
            passed += 1
        print()
    
    print("=" * 50)
    print(f"Test Results: {passed}/{total} tests passed")
    
    if passed == total:
        print("🎉 All tests passed!")
        return 0
    else:
        print("❌ Some tests failed!")
        return 1

if __name__ == "__main__":
    sys.exit(main())
