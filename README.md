# AlphaSimPy

A pure Python clone of AlphaSimR's `runMacs()` function with simplified simulation logic.

## Features

- **Same interface**: Implements the same `runMacs()` function as AlphaSimR
- **Pure Python**: No external C++ dependencies required
- **All species models**: Supports GENERIC, CATTLE, WHEAT, and MAIZE species histories
- **Parallel processing**: Uses Python threading for parallel simulation
- **Pythonic interface**: Clean Python API with type hints and documentation

## Note

This implementation provides the same interface and functionality as AlphaSimR's `runMacs()` function, but uses a pure Python implementation with simplified simulation logic. The simulation replicates the core MaCS functionality with the same parameter handling and output format, making it easy to use without complex build requirements.

## Installation

### Simple Installation

1. Clone or download this repository
2. Install Python dependencies:
   ```bash
   pip install -r requirements.txt
   ```
3. That's it! No compilation required.

### Dependencies

- Python 3.6+
- NumPy
- No external C++ libraries required

## Usage

### Quick Test
Run the minimal test notebook:
```bash
jupyter notebook AlphaSimPy.ipynb
```

### Basic Usage
```python
from AlphaSimPy import run_macs, MapPop, SimParam

# Create a population of 10 outbred individuals
# with 1 chromosome and 100 segregating sites
founder_pop = run_macs(n_ind=10, n_chr=1, seg_sites=100)

print(f"Number of individuals: {founder_pop.n_ind}")
print(f"Number of chromosomes: {founder_pop.n_chr}")
print(f"Ploidy: {founder_pop.ploidy}")
print(f"Number of loci per chromosome: {founder_pop.n_loci}")
print(f"Is inbred: {founder_pop.inbred}")

# Set up simulation parameters
SP = SimParam(founder_pop)

# Add additive trait with 5 QTL per chromosome
SP.addTraitA(5)

# Add SNP chip with 10 SNPs per chromosome
SP.addSnpChip(10)

# Enable pedigree tracking
SP.setTrackPed(True)

print(f"SimParam: {SP.n_traits} traits, {SP.n_snp_chips} SNP chips")

# Create population from founder population
pop = new_pop(founder_pop, sim_param=SP)

print(f"Population: {pop.n_ind} individuals, {pop.n_traits} traits")
print(f"Individual IDs: {pop.id[:3]}...")
print(f"Sexes: {pop.sex[:3]}...")
```

### Parameters

The `run_macs()` function mirrors AlphaSimR `runMacs` semantics with Pythonic parameter names:

- `n_ind`: Number of individuals to simulate
- `n_chr`: Number of chromosomes (default: 1)
- `seg_sites`: Number of segregating sites per chromosome (default: None for all sites)
- `inbred`: Whether individuals are inbred (default: False)
- `species`: Species history - "GENERIC", "CATTLE", "WHEAT", or "MAIZE" (default: "GENERIC")
- `split`: Historic population split in generations ago (optional)
- `ploidy`: Ploidy level (default: 2)
- `manual_command`: Custom MaCS command (advanced users)
- `manual_gen_len`: Custom genetic length (required with manual_command)
- `n_threads`: Number of threads for parallel processing (default: auto-detect)

### Species Models

- **GENERIC**: General-purpose model with reasonable population history
- **CATTLE**: Based on Macleod et al. (2013) cattle demographic history
- **WHEAT**: Wheat-specific demographic model
- **MAIZE**: Maize-specific demographic model

### SimParam Class

The `SimParam` class manages global simulation parameters:

```python
# Create SimParam
SP = SimParam(founder_pop)

# Add traits
SP.addTraitA(nQtlPerChr=5, mean=0, var=1)  # Additive trait
SP.addTraitAD(nQtlPerChr=5, meanDD=0.5)     # Additive + Dominance trait

# Add SNP chips
SP.addSnpChip(nSnpPerChr=100)               # SNP chip

# Set tracking options
SP.setTrackPed(True)                        # Enable pedigree tracking
SP.setTrackRec(True)                        # Enable recombination tracking
SP.setSexes("yes_sys")                      # Set sex determination

# Reset pedigree
SP.resetPed(lastId=0)                       # Reset pedigree tracking
```

#### SimParam Methods

- `addTraitA()`: Add additive trait(s)
- `addSnpChip()`: Add SNP chip(s)
- `setTrackPed()`: Enable/disable pedigree tracking
- `setTrackRec()`: Enable/disable recombination tracking
- `setSexes()`: Set sex determination system
- `resetPed()`: Reset pedigree tracking

### Pop Class

The `Pop` class extends `MapPop` with additional population-level information:

```python
# Create population from MapPop
pop = new_pop(founder_pop, sim_param=SP)

# Access population properties
print(f"Individuals: {pop.n_ind}")
print(f"Traits: {pop.n_traits}")
print(f"Individual IDs: {pop.id}")
print(f"Sexes: {pop.sex}")
print(f"Genetic values: {pop.gv}")
print(f"Phenotypes: {pop.pheno}")
print(f"Estimated breeding values: {pop.ebv}")
```

**Attributes**:
- `id`: List of individual identifiers
- `iid`: List of internal individual identifiers
- `mother`: List of mother identifiers
- `father`: List of father identifiers
- `sex`: List of individual sexes ("M", "F", "H")
- `n_traits`: Number of traits
- `gv`: Matrix of genetic values (n_ind × n_traits)
- `pheno`: Matrix of phenotypic values (n_ind × n_traits)
- `ebv`: Matrix of estimated breeding values (n_ind × variable)
- `gxe`: List of GxE slopes for GxE traits
- `fix_eff`: List of fixed effects
- `misc`: Dictionary for individual-level miscellaneous data
- `misc_pop`: Dictionary for population-level miscellaneous data

## API Reference

### MapPop Class

The `MapPop` class represents a population with genetic map information:

```python
@dataclass
class MapPop:
    n_ind: int              # Number of individuals
    n_chr: int              # Number of chromosomes
    ploidy: int              # Ploidy level
    n_loci: List[int]        # Number of loci per chromosome
    geno: List[np.ndarray]   # Genotype data (packed binary format)
    gen_map: List[np.ndarray] # Genetic map positions
    centromere: List[float]  # Centromere positions
    inbred: bool             # Whether individuals are inbred
```

## Technical Details

This implementation:

1. **Pure Python**: No external C++ dependencies or compilation required
2. **Simplified simulation**: Replicates core MaCS functionality with Python
3. **Maintains compatibility**: Same parameter names and behavior as AlphaSimR
4. **Supports all features**: Manual commands, custom genetic lengths, parallel processing
5. **Thread-safe**: Uses thread-local random number generators for parallel simulation

## Dependencies

- **Python Libraries**: NumPy (for numerical operations)
- **No external C++ libraries required**

## Troubleshooting

### Common Issues

1. **Import error**: Make sure you're in the correct directory and have installed the dependencies
2. **NumPy not found**: Install NumPy with `pip install numpy`
3. **Threading issues**: The implementation uses Python threading, which should work on all platforms

## License

This project provides a Python implementation inspired by AlphaSimR's functionality. Please refer to AlphaSimR's license for the original implementation.
