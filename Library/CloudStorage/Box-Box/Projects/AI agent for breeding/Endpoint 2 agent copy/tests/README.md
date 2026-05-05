# AlphaSimPy Tests

This directory contains unit tests for AlphaSimPy, translated from the AlphaSimR test suite.

## Test Files

- `test_phenotypes.py` - Tests for phenotype-related functions (asCategorical)
- `test_statistics.py` - Tests for statistical functions (addError)
- `test_hybrids.py` - Tests for hybrid crossing functions
- `test_misc.py` - Tests for miscellaneous functions
- `test_editGenome.py` - Tests for genome editing functions
- `test_crossing.py` - Tests for crossing functions
- `test_importData.py` - Tests for data import functions
- `test_selection.py` - Tests for selection functions
- `test_addTrait.py` - Tests for trait addition functions

## Running Tests

To run all tests:
```bash
pytest tests/
```

To run a specific test file:
```bash
pytest tests/test_phenotypes.py
```

## Notes

Some functions referenced in the tests may not yet be implemented in AlphaSimPy. These tests will fail until the corresponding functions are implemented. The tests are structured to match the original R testthat tests as closely as possible, with appropriate Python syntax conversions:

- R's `@` slot access becomes Python's `.` attribute access
- R's 1-based indexing becomes Python's 0-based indexing
- R's `expect_equal` becomes `np.testing.assert_array_equal` or `np.testing.assert_allclose`
- R's `expect_error` becomes `pytest.raises`
- R's `test_that` becomes Python's `def test_*` functions

## Functions That May Need Implementation

Based on the tests, the following functions may need to be implemented or exposed:

- `asCategorical` - Convert continuous values to categorical
- `addError` - Add error to genetic values (may be internal)
- `newMapPop` - Create MapPop from genetic map and haplotypes
- `quickHaplo` - Quick haplotype generation
- `editGenome` - Edit genome at specific loci
- `makeCross` - Make crosses from a cross plan
- `selectCross` - Select and cross individuals
- `selectOP` - Select open pollinated individuals
- `importHaplo` - Import haplotypes from external data
- `genParam` - Calculate genetic parameters
- `newMultiPop` - Create multi-population object
- `newEmptyPop` - Create empty population
- `pullSegSiteGeno` - Pull segregating site genotypes

