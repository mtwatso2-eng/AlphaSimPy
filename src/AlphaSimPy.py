import numpy as np
import random
import importlib
from typing import Optional, List, Union, Dict, Any
from dataclasses import dataclass
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
import threading

# Try to import C++ bindings, fall back to Python implementation if not available
try:
    try:
        import _alphasimpy_cpp
    except ImportError:
        # Package-style extension location
        _alphasimpy_cpp = importlib.import_module("alphasimpy._alphasimpy_cpp")
    _USE_CPP = True
except ImportError:
    _USE_CPP = False
    print("Warning: C++ bindings not available. Using Python implementation.")
    # Thread-local random number generator
    _thread_local = threading.local()

    def get_thread_local_random():
        """Get thread-local random number generator"""
        if not hasattr(_thread_local, 'rng'):
            _thread_local.rng = random.Random()
        return _thread_local.rng


@dataclass
class MapPop:
    """Python equivalent of AlphaSimR's MapPop class"""
    n_ind: int
    n_chr: int
    ploidy: int
    n_loci: List[int]
    geno: List[np.ndarray]
    gen_map: List[np.ndarray]
    centromere: List[float]
    inbred: bool
    gen_map_names: Optional[List[np.ndarray]] = None  # Marker names per chromosome
    
    def __post_init__(self):
        """Validate the MapPop object"""
        if self.n_chr != len(self.geno):
            raise ValueError("nChr != length(geno)")
        if self.n_chr != len(self.n_loci):
            raise ValueError("nChr != length(nLoci)")
        if self.n_chr != len(self.gen_map):
            raise ValueError("nChr != length(genMap)")
        if self.n_chr != len(self.centromere):
            raise ValueError("nChr != length(centromere)")
        
        for i in range(self.n_chr):
            expected_bins = self.n_loci[i] // 8 + (self.n_loci[i] % 8 > 0)
            if expected_bins != self.geno[i].shape[0]:
                raise ValueError(f"nLoci[{i}] != dim(geno[[{i}]][0])")
            if self.ploidy != self.geno[i].shape[1]:
                raise ValueError(f"ploidy != dim(geno[[{i}]][1])")
            if self.n_ind != self.geno[i].shape[2]:
                raise ValueError(f"nInd != dim(geno[[{i}]][2])")
            if self.n_loci[i] != len(self.gen_map[i]):
                raise ValueError(f"nLoci[{i}] != length(genMap[[{i}]])")


def get_num_threads() -> int:
    """Get the number of available threads"""
    if _USE_CPP:
        return _alphasimpy_cpp.get_num_threads()
    return mp.cpu_count()


def sample_int(n: int, N: int) -> np.ndarray:
    """Sample n integers without replacement from 0 to N-1"""
    if _USE_CPP:
        return _alphasimpy_cpp.sample_int(n, N)
    
    if n == 0:
        return np.array([], dtype=np.uint32)
    if n >= N:
        return np.arange(N, dtype=np.uint32)
    
    # Use numpy's random choice for efficiency
    return np.random.choice(N, size=n, replace=False).astype(np.uint32)


def _trans_mat(R: np.ndarray) -> np.ndarray:
    """
    Compute transformation matrix for correlated random normals.
    From AlphaSimR misc.R transMat - eigendecomposition for positive semi-definite check.
    """
    if not np.allclose(R, R.T):
        raise ValueError("Matrix is not symmetric")
    eig_vals, eig_vecs = np.linalg.eigh(R)
    eps = np.finfo(float).eps
    if np.min(eig_vals) < eps:
        import warnings
        warnings.warn("Matrix is not positive semi-definite, see ?transMat for details")
        eig_vals = np.maximum(eig_vals, 100 * eps)
        m = R.shape[0]
        tot_var = np.sum(eig_vals)
        eig_vals = eig_vals * m / tot_var
        new_R = eig_vecs @ np.diag(eig_vals) @ eig_vecs.T
        new_R = new_R / np.sqrt(np.outer(np.diag(new_R), np.diag(new_R)))
        eig_vals, eig_vecs = np.linalg.eigh(new_R)
    return (eig_vecs * np.sqrt(np.maximum(eig_vals, 0))) @ eig_vecs.T


def add_error(gv: np.ndarray, var_e: Union[np.ndarray, List[float]], reps: int = 1) -> np.ndarray:
    """
    Add residual error to genetic values.
    From AlphaSimR phenotypes.R addError.
    
    Parameters:
    -----------
    gv : numpy.ndarray
        Matrix of genetic values (n_ind x n_traits)
    var_e : numpy.ndarray or list
        Residual variances - vector (diagonal) or full covariance matrix
    reps : int, default=1
        Number of replications for phenotype (scales variance by 1/reps)
    
    Returns:
    --------
    numpy.ndarray
        Phenotype matrix (gv + error)
    """
    gv = np.asarray(gv)
    n_traits = gv.shape[1]
    n_ind = gv.shape[0]
    
    if isinstance(var_e, np.ndarray) and var_e.ndim == 2:
        if not np.allclose(var_e, var_e.T):
            raise ValueError("Covariance matrix must be symmetric")
        if var_e.shape[0] != n_traits:
            raise ValueError("Covariance matrix dimension must match n_traits")
        error = np.random.standard_normal((n_ind, n_traits)) @ _trans_mat(var_e)
    else:
        var_e = np.atleast_1d(var_e)
        if len(var_e) != n_traits:
            raise ValueError("var_e length must match n_traits")
        error_list = []
        for v in var_e:
            if np.isnan(v):
                error_list.append(np.full(n_ind, np.nan))
            else:
                error_list.append(np.random.normal(0, np.sqrt(float(v)), n_ind))
        error = np.column_stack(error_list)
    
    error = error / np.sqrt(reps)
    return gv + error


# Removed simulate_chromosome - now using C++ implementation via macs_wrapper


def run_macs(n_ind: int,
             n_chr: int = 1,
             seg_sites: Optional[Union[int, List[int]]] = None,
             inbred: bool = False,
             species: str = "GENERIC",
             split: Optional[int] = None,
             ploidy: int = 2,
             manual_command: Optional[str] = None,
             manual_gen_len: Optional[Union[float, List[float]]] = None,
             n_threads: Optional[int] = None) -> MapPop:
    """
    Create founder haplotypes using simplified MaCS-like simulation
    
    Parameters:
    -----------
    n_ind : int
        Number of individuals to simulate
    n_chr : int, default=1
        Number of chromosomes to simulate
    seg_sites : int or list of int, optional
        Number of segregating sites to keep per chromosome. If None, all sites are retained.
    inbred : bool, default=False
        Should founder individuals be inbred
    species : str, default="GENERIC"
        Species history to simulate. Options: "GENERIC", "CATTLE", "WHEAT", "MAIZE"
    split : int, optional
        Optional historic population split in terms of generations ago
    ploidy : int, default=2
        Ploidy level of organism
    manual_command : str, optional
        User provided MaCS options. For advanced users only.
    manual_gen_len : float or list of float, optional
        User provided genetic length. Must be supplied if using manual_command.
    n_threads : int, optional
        Number of threads for parallel simulation. If None, automatically detected.
    
    Returns:
    --------
    MapPop
        A MapPop object containing the simulated population
    """
    
    # Set default number of threads
    if n_threads is None:
        n_threads = get_num_threads()
    
    # Convert inputs to appropriate types
    n_ind = int(n_ind)
    n_chr = int(n_chr)
    ploidy = int(ploidy)
    
    if n_chr < n_threads:
        n_threads = n_chr
    
    # Generate random seeds for each chromosome
    seed = [str(random.randint(1, 100000000)) for _ in range(n_chr)]
    
    # Handle seg_sites parameter
    if seg_sites is None:
        seg_sites = [0] * n_chr
    elif isinstance(seg_sites, int):
        seg_sites = [int(seg_sites)] * n_chr
    else:
        seg_sites = [int(x) for x in seg_sites]
    
    # Calculate population size
    pop_size = n_ind if inbred else ploidy * n_ind
    
    # Build MaCS command
    if manual_command is not None:
        if manual_gen_len is None:
            raise ValueError("You must define manual_gen_len when using manual_command")
        command = f"{pop_size} {manual_command} -s "
        gen_len = manual_gen_len
    else:
        species = species.upper()
        
        if species == "GENERIC":
            gen_len = 1.0
            ne = 100
            species_params = "1E8 -t 1E-5 -r 4E-6"
            species_hist = "-eN 0.25 5.0 -eN 2.50 15.0 -eN 25.00 60.0 -eN 250.00 120.0 -eN 2500.00 1000.0"
            
        elif species == "CATTLE":
            # Macleod et al. (2013) https://doi.org/10.1093/molbev/mst125
            cattle_chr_sum = 2.8e9  # https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_002263795.3/
            cattle_chr_bp = cattle_chr_sum / 30
            rec_rate = 9.26e-9
            gen_len = rec_rate * cattle_chr_bp
            mut_rate = 9.4e-9
            ne = 90
            hist_ne = [120, 250, 350, 1000, 1500, 2000, 2500, 3500, 7000, 10000, 17000, 62000]
            hist_gen = [3, 6, 12, 18, 24, 154, 454, 654, 1754, 2354, 3354, 33154]
            species_params = f"{int(cattle_chr_bp)} -t {mut_rate * 4 * ne} -r {rec_rate * 4 * ne}"
            hist_ne = [x/ne for x in hist_ne]
            hist_gen = [x/(4*ne) for x in hist_gen]
            species_hist = ""
            for i in range(len(hist_ne)):
                species_hist += f" -eN {hist_gen[i]} {hist_ne[i]}"
                
        elif species == "WHEAT":
            gen_len = 1.43
            ne = 100
            species_params = "8E8 -t 4E-7 -r 3.6E-7"
            species_hist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 1.00 40 -eN 2.00 60 -eN 3.00 80 -eN 4.00 100 -eN 5.00 120 -eN 10.00 140 -eN 20.00 160 -eN 30.00 180 -eN 40.00 200 -eN 50.00 240 -eN 100.00 320 -eN 200.00 400 -eN 300.00 480 -eN 400.00 560 -eN 500.00 640"
            
        elif species == "MAIZE":
            gen_len = 2.0
            ne = 100
            species_params = "2E8 -t 5E-6 -r 4E-6"
            species_hist = "-eN 0.03 1 -eN 0.05 2 -eN 0.10 4 -eN 0.15 6 -eN 0.20 8 -eN 0.25 10 -eN 0.30 12 -eN 0.35 14 -eN 0.40 16 -eN 0.45 18 -eN 0.50 20 -eN 2.00 40 -eN 3.00 60 -eN 4.00 80 -eN 5.00 100"
            
        else:
            raise ValueError(f"No rules for species {species}")
        
        # Handle population split
        if split is None:
            split_i = ""
            split_j = ""
        else:
            if pop_size % 2 != 0:
                raise ValueError("Population size must be even when using split")
            split_i = f" -I 2 {pop_size//2} {pop_size//2}"
            split_j = f" -ej {split/(4*ne)+0.000001} 2 1"
        
        command = f"{pop_size} {species_params}{split_i} {species_hist}{split_j} -s "
    
    # Handle manual_gen_len override
    if manual_gen_len is not None:
        gen_len = manual_gen_len
    
    # Ensure gen_len is a list
    if isinstance(gen_len, (int, float)):
        gen_len = [float(gen_len)] * n_chr
    else:
        gen_len = [float(x) for x in gen_len]
    
    # Simulate chromosomes using C++ implementation
    if _USE_CPP:
        # Convert seg_sites to numpy array
        max_sites_array = np.array(seg_sites, dtype=np.uint32)
        
        # Call C++ MaCS function
        result = _alphasimpy_cpp.macs(
            command,
            max_sites_array,
            inbred,
            ploidy,
            n_threads,
            seed
        )
        
        geno_list = result['geno']
        gen_map_list = result['gen_map']
        n_loci_list = [len(gm) for gm in gen_map_list]
    else:
        # Fallback to Python implementation (removed - would need to be reimplemented)
        raise RuntimeError("C++ bindings required. Please build the extension module.")
    
    # Check if desired number of loci were obtained
    is_limited = [sites > 0 for sites in seg_sites]
    
    for i, (limited, expected, actual) in enumerate(zip(is_limited, seg_sites, n_loci_list)):
        if limited and actual != expected:
            raise ValueError(f"Simulation did not return enough segSites for chromosome {i+1}. "
                           f"Expected {expected}, got {actual}. Use seg_sites=None to return all sites.")
    
    # Convert numpy arrays back to lists for MapPop compatibility
    geno_list = [np.array(g) for g in geno_list]
    gen_map_list = [np.array(gm) for gm in gen_map_list]
    
    # Process genetic map
    processed_gen_map = []
    centromere = []
    
    for i in range(n_chr):
        chr_map = gen_map_list[i]
        # Scale by genetic length and normalize to start at 0
        scaled_map = gen_len[i] * (chr_map - chr_map[0])
        processed_gen_map.append(scaled_map)
        centromere.append(np.max(scaled_map) / 2)
    
    # Create MapPop object
    output = MapPop(
        n_ind=n_ind,
        n_chr=n_chr,
        ploidy=ploidy,
        n_loci=n_loci_list,
        geno=geno_list,
        gen_map=processed_gen_map,
        centromere=centromere,
        inbred=inbred
    )
    
    return output


def _pack_haplo(haplo: np.ndarray, ploidy: int, inbred: bool) -> np.ndarray:
    """
    Pack haplotype matrix into binary format. Matches AlphaSimR packHaplo.cpp logic.
    haplo: (n_haplo, n_loci) matrix of 0/1
    Returns: (n_bins, ploidy, n_ind) uint8 array
    """
    haplo = np.asarray(haplo, dtype=np.uint8)
    n_hap, n_loci = haplo.shape
    if inbred:
        n_ind = n_hap
    else:
        if n_hap % ploidy != 0:
            raise ValueError("Number of rows not a factor of ploidy")
        n_ind = n_hap // ploidy
    n_bins = n_loci // 8 + (1 if n_loci % 8 else 0)
    output = np.zeros((n_bins, ploidy, n_ind), dtype=np.uint8)
    for i in range(n_hap):
        locus = 0
        for j in range(n_bins):
            byte_val = 0
            for k in range(8):
                if locus < n_loci:
                    byte_val |= (int(haplo[i, locus]) << k)
                    locus += 1
            if inbred:
                for chr_idx in range(ploidy):
                    output[j, chr_idx, i] = byte_val
            else:
                chr_idx = i % ploidy
                ind_idx = i // ploidy
                output[j, chr_idx, ind_idx] = byte_val
    return output


def new_map_pop(gen_map: List[np.ndarray], haplotypes: List[np.ndarray],
                inbred: bool = False, ploidy: int = 2) -> MapPop:
    """
    Create MapPop from genetic map and haplotypes. From AlphaSimR founderPop.R newMapPop.
    """
    ploidy = int(ploidy)
    if len(gen_map) != len(haplotypes):
        raise ValueError("genMap and haplotypes must have same length")
    n_row = [h.shape[0] for h in haplotypes]
    if len(set(n_row)) > 1:
        raise ValueError("Number of rows must be equal in haplotypes")
    n_row = n_row[0]
    if inbred:
        n_ind = n_row
    else:
        if n_row % ploidy != 0:
            raise ValueError("Number of haplotypes must be divisible by ploidy")
        n_ind = n_row // ploidy
    seg_sites = [len(g) for g in gen_map]
    n_col = [h.shape[1] for h in haplotypes]
    if not all(nc == seg_sites[i] for i, nc in enumerate(n_col)):
        raise ValueError("Number of segregating sites in haplotypes and genMap don't match")
    geno_list = []
    for chr_idx in range(len(gen_map)):
        packed = _pack_haplo(haplotypes[chr_idx], ploidy, inbred)
        geno_list.append(packed)
    centromere = [np.max(g) / 2 for g in gen_map]
    return MapPop(
        n_ind=n_ind, n_chr=len(gen_map), ploidy=ploidy,
        n_loci=seg_sites, geno=geno_list, gen_map=gen_map,
        centromere=centromere, inbred=inbred
    )


def quick_haplo(n_ind: int, n_chr: int, seg_sites: Union[int, List[int]],
                gen_len: Union[float, List[float]] = 1, ploidy: int = 2,
                inbred: bool = False) -> MapPop:
    """
    Quick founder haplotype simulation - random 0/1 sampling.
    From AlphaSimR founderPop.R quickHaplo.
    """
    ploidy = int(ploidy)
    n_ind = int(n_ind)
    n_chr = int(n_chr)
    seg_sites = [int(seg_sites)] * n_chr if np.isscalar(seg_sites) else [int(s) for s in seg_sites]
    gen_len = [float(gen_len)] * n_chr if np.isscalar(gen_len) else [float(g) for g in gen_len]
    n_bins = [s // 8 + (1 if s % 8 else 0) for s in seg_sites]
    centromere = [g / 2 for g in gen_len]
    geno_list = []
    gen_map_list = []
    for i in range(n_chr):
        gen_map_list.append(np.linspace(0, gen_len[i], seg_sites[i]))
        raw_geno = np.random.randint(0, 256, size=(n_ind * ploidy * n_bins[i]), dtype=np.uint8)
        geno_arr = raw_geno.reshape(n_bins[i], ploidy, n_ind)
        if inbred and ploidy > 1:
            for j in range(1, ploidy):
                geno_arr[:, j, :] = geno_arr[:, 0, :]
        geno_list.append(geno_arr)
    return MapPop(
        n_ind=n_ind, n_chr=n_chr, ploidy=ploidy,
        n_loci=seg_sites, geno=geno_list, gen_map=gen_map_list,
        centromere=centromere, inbred=inbred
    )


def import_gen_map(gen_map) -> tuple:
    """
    Format genetic map from DataFrame to (positions_dict, names_dict).
    From AlphaSimR importData.R importGenMap.
    gen_map: DataFrame with columns markerName, chromosome, position
    """
    try:
        import pandas as pd
    except ImportError:
        raise ImportError("pandas required for importGenMap")
    if isinstance(gen_map, pd.DataFrame):
        marker_name = gen_map.iloc[:, 0].astype(str).values
        chromosome = gen_map.iloc[:, 1].astype(str).values
        position = gen_map.iloc[:, 2].astype(float).values
    else:
        arr = np.asarray(gen_map)
        marker_name = arr[:, 0].astype(str)
        chromosome = arr[:, 1].astype(str)
        position = arr[:, 2].astype(float)
    unique_chr = np.unique(chromosome)
    gen_map_pos = {}
    gen_map_names = {}
    for ch in unique_chr:
        take = chromosome == ch
        pos = position[take].copy()
        names = marker_name[take].copy()
        order = np.argsort(pos)
        pos = pos[order]
        names = names[order]
        pos = pos - pos[0]
        gen_map_pos[str(ch)] = pos
        gen_map_names[str(ch)] = names
    return gen_map_pos, gen_map_names


def import_haplo(haplo, gen_map, ploidy: int = 2, ped=None) -> MapPop:
    """
    Import haplotypes from matrix format. From AlphaSimR importData.R importHaplo.
    """
    gen_map_pos, gen_map_names = import_gen_map(gen_map)
    haplo = np.asarray(haplo)
    if haplo.ndim == 1:
        haplo = haplo.reshape(1, -1)
    haplo = haplo.astype(np.uint8)
    if not np.all((haplo == 0) | (haplo == 1)):
        raise ValueError("Haplo must contain only 0 and 1")
    if isinstance(gen_map, np.ndarray):
        all_marker_names = gen_map[:, 0].astype(str)
    else:
        try:
            import pandas as pd
            all_marker_names = gen_map.iloc[:, 0].astype(str).values
        except Exception:
            all_marker_names = np.asarray(gen_map)[:, 0].astype(str)
    chr_names = sorted(gen_map_pos.keys(), key=lambda x: (int(x) if str(x).isdigit() else x))
    haplotypes = []
    gen_map_arrays = []
    names_arrays = []
    for ch in chr_names:
        map_markers = gen_map_names[ch]
        pos_arr = gen_map_pos[ch]
        take = []
        valid_pos = []
        valid_names = []
        for i, m in enumerate(map_markers):
            idx = np.where(all_marker_names == m)[0]
            if len(idx) > 0:
                take.append(idx[0])
                valid_pos.append(pos_arr[i])
                valid_names.append(m)
        if len(take) < 1:
            raise ValueError("No matching markers found for chromosome")
        haplotypes.append(haplo[:, take])
        gen_map_arrays.append(np.asarray(valid_pos))
        names_arrays.append(np.array(valid_names))
    founder = new_map_pop(gen_map_arrays, haplotypes, ploidy=ploidy)
    founder.gen_map_names = names_arrays
    return founder


def edit_genome(pop, ind, chr, seg_sites, allele, sim_param=None) -> 'Pop':
    """
    Edit genome at specific locus. From AlphaSimR misc.R editGenome.
    Sets the specified allele at the given segregating site for the given individual(s).
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    ind = np.unique(np.atleast_1d(ind).astype(int))
    if all(1 <= i <= pop.n_ind for i in ind):
        ind = ind - 1
    elif all(0 <= i < pop.n_ind for i in ind):
        ind = ind
    else:
        raise ValueError("Invalid individual index")
    chr = np.atleast_1d(chr).astype(int)
    seg_sites = np.atleast_1d(seg_sites).astype(int)
    if len(chr) != len(seg_sites):
        raise ValueError("chr and segSites must have same length")
    allele = np.atleast_1d(allele).astype(int)
    if not np.all((allele == 0) | (allele == 1)):
        raise ValueError("allele must be 0 or 1")
    if len(allele) == 1:
        allele = np.repeat(allele, len(seg_sites))
    if len(allele) != len(seg_sites):
        raise ValueError("allele length must match segSites")
    allele = allele.astype(np.uint8)
    for sel_chr in np.unique(chr):
        sel = np.where(chr == sel_chr)[0]
        for i in sel:
            byte_idx = (seg_sites[i] - 1) // 8
            bit_idx = (seg_sites[i] - 1) % 8
            for sel_ind in ind:
                idx = sel_ind
                for j in range(pop.ploidy):
                    old_byte = pop.geno[sel_chr - 1][byte_idx, j, idx]
                    bits = [(old_byte >> k) & 1 for k in range(8)]
                    bits[bit_idx] = allele[i]
                    new_byte = sum(b << k for k, b in enumerate(bits))
                    pop.geno[sel_chr - 1][byte_idx, j, idx] = new_byte
    gv_new = _get_gv_index(pop, sim_param)
    pop.gv = gv_new
    pop.pheno = gv_new.copy()
    return pop


def new_empty_pop(ploidy: int = 2, sim_param=None) -> 'Pop':
    """Create empty population. From AlphaSimR Class-Pop.R newEmptyPop."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    n_loci = sim_param.founder_pop.n_loci
    geno = []
    for i in range(sim_param.n_chr):
        n_bins = n_loci[i] // 8 + (1 if n_loci[i] % 8 else 0)
        geno.append(np.zeros((n_bins, ploidy, 0), dtype=np.uint8))
    return Pop(
        n_ind=0, n_chr=sim_param.n_chr, ploidy=ploidy, n_loci=n_loci,
        geno=geno, gen_map=sim_param.founder_pop.gen_map,
        centromere=sim_param.founder_pop.centromere, inbred=False,
        id=[], iid=[], mother=[], father=[], sex=[],
        n_traits=sim_param.n_traits,
        gv=np.empty((0, sim_param.n_traits)), pheno=np.empty((0, sim_param.n_traits)),
        ebv=np.empty((0, 0)), gxe=[None] * sim_param.n_traits,
        fix_eff=[], misc={}, misc_pop={}
    )


@dataclass
class MultiPop:
    """Container for multiple Pop objects. From AlphaSimR MultiPop."""
    pops: List['Pop']
    
    def __len__(self):
        return len(self.pops)
    
    def __getitem__(self, i):
        if isinstance(i, slice):
            return MultiPop(pops=self.pops[i])
        return self.pops[i]


def new_multi_pop(*pops) -> MultiPop:
    """Create MultiPop from one or more Pop objects. From AlphaSimR newMultiPop."""
    input_list = list(pops)
    for p in input_list:
        if not isinstance(p, (Pop, MultiPop)):
            raise ValueError("All arguments must be Pop or MultiPop")
    flat = []
    for p in input_list:
        if isinstance(p, MultiPop):
            flat.extend(p.pops)
        else:
            flat.append(p)
    return MultiPop(pops=flat)


def _pull_geno_from_packed(geno_list: List[np.ndarray], n_loci: List[int],
                           ploidy: int, loci_per_chr: List[np.ndarray]) -> np.ndarray:
    """Extract genotype matrix from packed format for specified loci."""
    all_geno = []
    for chr_idx, loci in enumerate(loci_per_chr):
        if len(loci) == 0:
            continue
        geno = geno_list[chr_idx]
        n_ind = geno.shape[2]
        for loc in loci:
            byte_idx = loc // 8
            bit_idx = loc % 8
            dosage = np.zeros(n_ind, dtype=int)
            for p in range(ploidy):
                byte_val = geno[byte_idx, p, :]
                dosage += (byte_val >> bit_idx) & 1
            all_geno.append(dosage)
    return np.column_stack(all_geno) if all_geno else np.empty((geno_list[0].shape[2], 0))


def pull_seg_site_geno(pop, chr=None, as_raw: bool = False,
                       sim_param=None) -> np.ndarray:
    """Pull segregating site genotypes. From AlphaSimR pullGeno.R pullSegSiteGeno."""
    if hasattr(pop, 'geno'):
        n_loci = list(pop.n_loci)
        geno = pop.geno
        ploidy = pop.ploidy
    else:
        if sim_param is None:
            raise ValueError("simParam must be provided for MapPop")
        n_loci = list(sim_param._seg_sites)
        geno = pop.geno
        ploidy = pop.ploidy
    if chr is not None:
        chr_idx = np.atleast_1d(chr) - 1
        geno = [geno[i] for i in chr_idx if 0 <= i < len(geno)]
        n_loci = [n_loci[i] for i in chr_idx if 0 <= i < len(n_loci)]
    loci_per_chr = [np.arange(n) for n in n_loci]
    result = _pull_geno_from_packed(geno, n_loci, ploidy, loci_per_chr)
    return result


# Trait classes
@dataclass
class TraitA:
    """Additive trait class"""
    qtl_loci: List[int]
    add_eff: np.ndarray
    intercept: float
    name: str
    
    @property
    def addEff(self):
        return self.add_eff
    
    @property
    def nLoci(self):
        return len(self.qtl_loci)


@dataclass
class TraitAG:
    """Additive + GxE trait class"""
    qtl_loci: List[int]
    add_eff: np.ndarray
    intercept: float
    name: str
    gxe_eff: np.ndarray
    gxe_int: float
    env_var: float


@dataclass
class SNP:
    """SNP marker class"""
    loci: List[int]
    name: str


class SimParam:
    """
    Simulation parameters container for AlphaSimPy
    
    This class manages global simulation parameters including traits, SNP chips,
    pedigree tracking, and other simulation settings.
    """
    
    def __init__(self, founder_pop: MapPop = None, **kwargs):
        """
        Initialize SimParam with a founder population
        
        Parameters:
        -----------
        founder_pop : MapPop
            The founder population used for variance scaling
        """
        founder_pop = founder_pop or kwargs.get('founder_pop')
        if founder_pop is None:
            raise ValueError("founder_pop must be provided")
        # Public attributes
        self.n_threads = get_num_threads()
        self.v = 2.6  # Kosambi mapping function
        self.p = 0    # Single pathway gamma model
        self.quad_prob = 0  # No quadrivalent pairing
        self.snp_chips = []
        self.invalid_qtl = [[] for _ in range(founder_pop.n_chr)]  # All eligible
        self.invalid_snp = [[] for _ in range(founder_pop.n_chr)]  # All eligible
        self.founder_pop = founder_pop
        self.finalize_pop = lambda pop: pop  # Default function
        self.allow_empty_pop = False  # Empty populations trigger an error
        
        # Private attributes
        self._restr_sites = True
        self._traits = []
        self._seg_sites = founder_pop.n_loci
        self._sexes = "no"
        self._female_map = [gm - gm[0] for gm in founder_pop.gen_map]  # Set position 1 to 0
        self._male_map = None
        self._sep_map = False
        self._female_centromere = founder_pop.centromere
        self._male_centromere = None
        self._last_id = 0
        self._is_track_ped = False
        self._pedigree = np.empty((0, 3), dtype=int)
        self._is_track_rec = False
        self._rec_hist = []
        self._var_a = []
        self._var_g = []
        self._var_e = []
        self._version = "1.0.0"  # AlphaSimPy version
        self._last_haplo = 0
        self._has_hap = []
        self._hap = []
        self._is_founder = []
        self._active_qtl = None
        self._qtl_index = []
    
    @property
    def n_chr(self) -> int:
        """Number of chromosomes"""
        return self.founder_pop.n_chr
    
    @property
    def n_traits(self) -> int:
        """Number of traits"""
        return len(self._traits)
    
    @property
    def traits(self):
        """List of traits (AlphaSimR compatibility)"""
        return self._traits
    
    @property
    def n_snp_chips(self) -> int:
        """Number of SNP chips"""
        return len(self.snp_chips)
    
    def set_track_ped(self, is_track_ped: bool, force: bool = False):
        """
        Set pedigree tracking for the simulation
        
        Parameters:
        -----------
        is_track_ped : bool
            Should pedigree tracking be on
        force : bool, default=False
            Should the check for a running simulation be ignored
        """
        if not isinstance(is_track_ped, bool):
            raise ValueError("is_track_ped must be a boolean")
        
        if not force:
            self._is_running()
        
        self._is_track_ped = is_track_ped
        if not is_track_ped:
            self._is_track_rec = False
        
        return self
    
    def set_track_rec(self, is_track_rec: bool, force: bool = False):
        """
        Set recombination tracking for the simulation
        
        Parameters:
        -----------
        is_track_rec : bool
            Should recombination tracking be on
        force : bool, default=False
            Should the check for a running simulation be ignored
        """
        if not isinstance(is_track_rec, bool):
            raise ValueError("is_track_rec must be a boolean")
        
        if not force:
            self._is_running()
        
        self._is_track_rec = is_track_rec
        if is_track_rec:
            self._is_track_ped = True
        
        return self
    
    def reset_ped(self, last_id: int = 0):
        """
        Reset the internal lastId, pedigree and recombination tracking
        
        Parameters:
        -----------
        last_id : int, default=0
            Last ID to include in pedigree
        """
        self._last_id = last_id
        if last_id > 0:
            self._pedigree = self._pedigree[:last_id]
        else:
            self._pedigree = np.empty((0, 3), dtype=int)
        
        if self._is_track_rec:
            self._rec_hist = self._rec_hist[:last_id]
        
        return self
    
    def set_sexes(self, sexes: str, force: bool = False):
        """
        Set sex determination system
        
        Parameters:
        -----------
        sexes : str
            Sex determination system. Options: "no", "yes_sys", "yes_rand"
        force : bool, default=False
            Should the check for a running simulation be ignored
        """
        if not force:
            self._is_running()
        
        if sexes not in ["no", "yes_sys", "yes_rand"]:
            raise ValueError("sexes must be one of: 'no', 'yes_sys', 'yes_rand'")
        
        self._sexes = sexes
        return self
    
    def add_trait_a(self, n_qtl_per_chr: Union[int, List[int]], 
                   mean: Union[float, List[float]] = 0,
                   var: Union[float, List[float]] = 1,
                   cor_a: Optional[np.ndarray] = None,
                   gamma: Union[bool, List[bool]] = False,
                   shape: Union[float, List[float]] = 1,
                   force: bool = False,
                   name: Optional[Union[str, List[str]]] = None):
        """
        Add additive trait(s) to the simulation
        
        Parameters:
        -----------
        n_qtl_per_chr : int or list of int
            Number of QTLs per chromosome
        mean : float or list of float, default=0
            Desired mean genetic values
        var : float or list of float, default=1
            Desired genetic variances
        cor_a : numpy.ndarray, optional
            Correlation matrix between additive effects
        gamma : bool or list of bool, default=False
            Should gamma distribution be used instead of normal
        shape : float or list of float, default=1
            Shape parameter for gamma distribution
        force : bool, default=False
            Should the check for a running simulation be ignored
        name : str or list of str, optional
            Optional names for traits
        """
        if not force:
            self._is_running()
        
        # Handle single values
        if isinstance(n_qtl_per_chr, int):
            n_qtl_per_chr = [n_qtl_per_chr] * self.n_chr
        
        if isinstance(mean, (int, float)):
            mean = [mean]
        if isinstance(var, (int, float)):
            var = [var]
        if isinstance(gamma, bool):
            gamma = [gamma] * len(mean)
        if isinstance(shape, (int, float)):
            shape = [shape] * len(mean)
        
        n_traits = len(mean)
        
        # Set up correlation matrix
        if cor_a is None:
            cor_a = np.eye(n_traits)
        
        # Set up trait names
        if name is None:
            name = [f"Trait{i+self.n_traits+1}" for i in range(n_traits)]
        elif isinstance(name, str):
            name = [name]
        
        # Validate inputs
        if len(mean) != len(var):
            raise ValueError("mean and var must have the same length")
        if not np.allclose(cor_a, cor_a.T):
            raise ValueError("cor_a must be symmetric")
        if cor_a.shape[0] != n_traits:
            raise ValueError("cor_a must be square with dimension equal to number of traits")
        if len(name) != n_traits:
            raise ValueError("name must have same length as mean")
        
        # Pick QTL loci
        qtl_loci = self._pick_loci(n_qtl_per_chr)
        
        # Flatten QTL loci from list of lists to flat list
        flat_qtl_loci = []
        for chr_loci in qtl_loci:
            flat_qtl_loci.extend(chr_loci)
        
        # Sample additive effects
        add_eff = self._samp_add_eff(qtl_loci, n_traits, cor_a, gamma, shape)
        
        # Create traits
        for i in range(n_traits):
            trait = TraitA(
                qtl_loci=flat_qtl_loci,
                add_eff=add_eff[:, i],
                intercept=0,
                name=name[i]
            )
            
            # Calculate genetic parameters
            tmp = self._calc_gen_param(trait, self.founder_pop)
            bv_var = self._pop_var(tmp['bv'])
            if np.isscalar(bv_var):
                scale = np.sqrt(var[i]) / np.sqrt(bv_var)
            else:
                scale = np.sqrt(var[i]) / np.sqrt(bv_var[0])
            trait.add_eff = trait.add_eff * scale
            trait.intercept = mean[i] - np.mean(tmp['gv'] * scale)
            
            self._add_trait(trait, var[i], var[i])
        
        return self
    
    def restr_seg_sites(self, min_qtl_per_chr: Optional[Union[int, List[int]]] = None,
                       min_snp_per_chr: Optional[Union[int, List[int]]] = None,
                       exclude_qtl: Optional[List[str]] = None,
                       exclude_snp: Optional[List[str]] = None,
                       overlap: bool = False,
                       min_snp_freq: Optional[float] = None):
        """
        Set restrictions on which segregating sites can serve as SNP and/or QTL
        
        Parameters:
        -----------
        min_qtl_per_chr : int or list of int, optional
            Minimum number of segregating sites for QTLs
        min_snp_per_chr : int or list of int, optional
            Minimum number of segregating sites for SNPs
        exclude_qtl : list of str, optional
            Segregating site names to exclude from QTL consideration
        exclude_snp : list of str, optional
            Segregating site names to exclude from SNP consideration
        overlap : bool, default=False
            Should SNP and QTL sites be allowed to overlap
        min_snp_freq : float, optional
            Minimum allowable frequency for SNP loci
        """
        # Handle single values
        if isinstance(min_qtl_per_chr, int):
            min_qtl_per_chr = [min_qtl_per_chr] * self.n_chr
        if isinstance(min_snp_per_chr, int):
            min_snp_per_chr = [min_snp_per_chr] * self.n_chr
        
        # Set overlap behavior
        self._restr_sites = not overlap
        
        # Pre-allocate sites if specified
        if min_qtl_per_chr is not None:
            for chr_idx in range(self.n_chr):
                n_qtl = min_qtl_per_chr[chr_idx]
                if n_qtl > 0:
                    # Randomly select QTL sites
                    available = [i for i in range(self._seg_sites[chr_idx]) 
                               if i not in self.invalid_qtl[chr_idx]]
                    if len(available) >= n_qtl:
                        selected = np.random.choice(available, size=n_qtl, replace=False)
                        self.invalid_qtl[chr_idx].extend(selected)
                        if not overlap:
                            self.invalid_snp[chr_idx].extend(selected)
        
        if min_snp_per_chr is not None:
            for chr_idx in range(self.n_chr):
                n_snp = min_snp_per_chr[chr_idx]
                if n_snp > 0:
                    # Randomly select SNP sites
                    available = [i for i in range(self._seg_sites[chr_idx]) 
                               if i not in self.invalid_snp[chr_idx]]
                    if len(available) >= n_snp:
                        selected = np.random.choice(available, size=n_snp, replace=False)
                        self.invalid_snp[chr_idx].extend(selected)
                        if not overlap:
                            self.invalid_qtl[chr_idx].extend(selected)
        
        return self
    
    def add_trait_ag(self, n_qtl_per_chr: Union[int, List[int]], 
                    mean: Union[float, List[float]] = 0,
                    var: Union[float, List[float]] = 1,
                    var_gxe: Union[float, List[float]] = 1e-6,
                    var_env: Union[float, List[float]] = 0,
                    cor_a: Optional[np.ndarray] = None,
                    cor_gxe: Optional[np.ndarray] = None,
                    gamma: Union[bool, List[bool]] = False,
                    shape: Union[float, List[float]] = 1,
                    force: bool = False,
                    name: Optional[Union[str, List[str]]] = None):
        """
        Add additive + GxE trait(s) to the simulation
        
        Parameters:
        -----------
        n_qtl_per_chr : int or list of int
            Number of QTLs per chromosome
        mean : float or list of float, default=0
            Desired mean genetic values
        var : float or list of float, default=1
            Desired genetic variances
        var_gxe : float or list of float, default=1e-6
            Total genotype-by-environment variances
        var_env : float or list of float, default=0
            Environmental variances
        cor_a : numpy.ndarray, optional
            Correlation matrix between additive effects
        cor_gxe : numpy.ndarray, optional
            Correlation matrix between GxE effects
        gamma : bool or list of bool, default=False
            Should gamma distribution be used instead of normal
        shape : float or list of float, default=1
            Shape parameter for gamma distribution
        force : bool, default=False
            Should the check for a running simulation be ignored
        name : str or list of str, optional
            Optional names for traits
        """
        if not force:
            self._is_running()
        
        # Handle single values
        if isinstance(n_qtl_per_chr, int):
            n_qtl_per_chr = [n_qtl_per_chr] * self.n_chr
        
        if isinstance(mean, (int, float)):
            mean = [mean]
        if isinstance(var, (int, float)):
            var = [var]
        if isinstance(var_gxe, (int, float)):
            var_gxe = [var_gxe]
        if isinstance(var_env, (int, float)):
            var_env = [var_env]
        if isinstance(gamma, bool):
            gamma = [gamma] * len(mean)
        if isinstance(shape, (int, float)):
            shape = [shape] * len(mean)
        
        n_traits = len(mean)
        
        # Set up correlation matrices
        if cor_a is None:
            cor_a = np.eye(n_traits)
        if cor_gxe is None:
            cor_gxe = np.eye(n_traits)
        
        # Set up trait names
        if name is None:
            name = [f"Trait{i+self.n_traits+1}" for i in range(n_traits)]
        elif isinstance(name, str):
            name = [name]
        
        # Validate inputs
        if len(mean) != len(var):
            raise ValueError("mean and var must have the same length")
        if len(var_gxe) != n_traits:
            raise ValueError("var_gxe must have same length as mean")
        if len(var_env) != n_traits:
            raise ValueError("var_env must have same length as mean")
        if not np.allclose(cor_a, cor_a.T):
            raise ValueError("cor_a must be symmetric")
        if not np.allclose(cor_gxe, cor_gxe.T):
            raise ValueError("cor_gxe must be symmetric")
        if cor_a.shape[0] != n_traits:
            raise ValueError("cor_a must be square with dimension equal to number of traits")
        if cor_gxe.shape[0] != n_traits:
            raise ValueError("cor_gxe must be square with dimension equal to number of traits")
        if len(name) != n_traits:
            raise ValueError("name must have same length as mean")
        
        # Pick QTL loci
        qtl_loci = self._pick_loci(n_qtl_per_chr)
        
        # Flatten QTL loci from list of lists to flat list
        flat_qtl_loci = []
        for chr_loci in qtl_loci:
            flat_qtl_loci.extend(chr_loci)
        
        # Sample additive effects
        add_eff = self._samp_add_eff(qtl_loci, n_traits, cor_a, gamma, shape)
        
        # Sample GxE effects
        gxe_eff = self._samp_add_eff(qtl_loci, n_traits, cor_gxe, gamma=[False] * n_traits, shape=shape)
        
        # Create traits
        for i in range(n_traits):
            # Create base additive trait
            trait = TraitA(
                qtl_loci=flat_qtl_loci,
                add_eff=add_eff[:, i],
                intercept=0,
                name=name[i]
            )
            
            # Calculate genetic parameters for additive effects
            tmp = self._calc_gen_param(trait, self.founder_pop)
            bv_var = self._pop_var(tmp['bv'])
            if np.isscalar(bv_var):
                scale = np.sqrt(var[i]) / np.sqrt(bv_var)
            else:
                scale = np.sqrt(var[i]) / np.sqrt(bv_var[0])
            trait.add_eff = trait.add_eff * scale
            trait.intercept = mean[i] - np.mean(tmp['gv'] * scale)
            
            # Calculate genetic parameters for GxE effects
            trait_gxe = TraitA(
                qtl_loci=flat_qtl_loci,
                add_eff=gxe_eff[:, i],
                intercept=0,
                name=name[i]
            )
            tmp_gxe = self._calc_gen_param(trait_gxe, self.founder_pop)
            
            # Create TraitAG
            if var_env[i] == 0:
                scale_gxe = np.sqrt(var_gxe[i]) / np.sqrt(self._pop_var(tmp_gxe['gv']))
                trait_ag = TraitAG(
                    qtl_loci=flat_qtl_loci,
                    add_eff=trait.add_eff,
                    intercept=trait.intercept,
                    name=name[i],
                    gxe_eff=gxe_eff[:, i] * scale_gxe,
                    gxe_int=0 - np.mean(tmp_gxe['gv'] * scale_gxe),
                    env_var=1
                )
            else:
                scale_gxe = np.sqrt(var_gxe[i] / var_env[i]) / np.sqrt(self._pop_var(tmp_gxe['gv']))
                trait_ag = TraitAG(
                    qtl_loci=flat_qtl_loci,
                    add_eff=trait.add_eff,
                    intercept=trait.intercept,
                    name=name[i],
                    gxe_eff=gxe_eff[:, i] * scale_gxe,
                    gxe_int=1 - np.mean(tmp_gxe['gv'] * scale_gxe),
                    env_var=var_env[i]
                )
            
            self._add_trait(trait_ag, var[i], var[i])
        
        return self
    
    def add_snp_chip(self, n_snp_per_chr: Union[int, List[int]], 
                    name: Optional[str] = None):
        """
        Add SNP chip to the simulation
        
        Parameters:
        -----------
        n_snp_per_chr : int or list of int
            Number of SNPs per chromosome
        name : str, optional
            Name for the SNP chip
        """
        if isinstance(n_snp_per_chr, int):
            n_snp_per_chr = [n_snp_per_chr] * self.n_chr
        
        # Pick SNP loci
        snp_loci = self._pick_loci(n_snp_per_chr, is_qtl=False)
        
        if name is None:
            name = f"SNPChip{self.n_snp_chips + 1}"
        
        snp_chip = SNP(loci=snp_loci, name=name)
        self.snp_chips.append(snp_chip)
        
        return self
    
    def import_trait(self, marker_names, add_eff, dom_eff=None, intercept=None,
                    name=None, var_e=None, force=False):
        """
        Import trait from marker effects. From AlphaSimR SimParam importTrait.
        """
        add_eff = np.atleast_2d(np.asarray(add_eff))
        if add_eff.shape[0] == 1:
            add_eff = add_eff.T
        marker_names = np.atleast_1d(marker_names)
        if len(marker_names) != add_eff.shape[0]:
            raise ValueError("markerNames length must match addEff rows")
        n_traits = add_eff.shape[1]
        intercept = np.zeros(n_traits) if intercept is None else np.atleast_1d(intercept)
        if len(intercept) != n_traits:
            raise ValueError("intercept length must match n_traits")
        name = name or [f"Trait{i+self.n_traits+1}" for i in range(n_traits)]
        if not isinstance(name, (list, tuple, np.ndarray)):
            name = [name]
        if len(name) != n_traits:
            raise ValueError("name length must match n_traits")
        gen_map_names = getattr(self.founder_pop, 'gen_map_names', None)
        if gen_map_names is None:
            raise ValueError("Founder population must have gen_map_names for importTrait (use importHaplo)")
        flat_loci = []
        flat_eff = []
        offset = 0
        for chr_idx in range(self.n_chr):
            chr_names = gen_map_names[chr_idx]
            for loc_idx, m in enumerate(chr_names):
                if m in marker_names:
                    pos = np.where(marker_names == m)[0][0]
                    flat_loci.append(offset + loc_idx)
                    flat_eff.append(add_eff[pos, :])
            offset += self.founder_pop.n_loci[chr_idx]
        if len(flat_loci) == 0:
            raise ValueError("No marker names matched the genetic map")
        flat_eff = np.array(flat_eff)
        for trt_idx in range(n_traits):
            trait = TraitA(qtl_loci=flat_loci, add_eff=flat_eff[:, trt_idx],
                          intercept=intercept[trt_idx], name=name[trt_idx])
            var_a = np.var(flat_eff[:, trt_idx])
            self._add_trait(trait, var_a, var_a)
        return self
    
    def importTrait(self, markerNames, addEff, domEff=None, intercept=None, name=None, varE=None, force=False):
        """Convenience method with AlphaSimR-style parameter names"""
        return self.import_trait(markerNames, addEff, domEff, intercept, name, varE, force)
    
    def _is_running(self):
        """Check if simulation is running (placeholder)"""
        # In a full implementation, this would check if populations are being created
        pass
    
    def _pick_loci(self, n_loci_per_chr: List[int], is_qtl: bool = True) -> List[List[int]]:
        """
        Pick loci for QTL or SNP
        
        Parameters:
        -----------
        n_loci_per_chr : list of int
            Number of loci per chromosome
        is_qtl : bool, default=True
            Whether picking QTL (True) or SNP (False) loci
        
        Returns:
        --------
        list of list of int
            Selected loci for each chromosome
        """
        invalid_loci = self.invalid_qtl if is_qtl else self.invalid_snp
        loci = []
        
        for chr_idx in range(self.n_chr):
            n_loci = n_loci_per_chr[chr_idx]
            n_total = self._seg_sites[chr_idx]
            
            # Get available loci
            available = [i for i in range(n_total) if i not in invalid_loci[chr_idx]]
            
            if len(available) < n_loci:
                raise ValueError(f"Not enough available loci on chromosome {chr_idx + 1}")
            
            # Randomly sample loci
            selected = np.random.choice(available, size=n_loci, replace=False)
            loci.append(sorted(selected.tolist()))
            
            # Mark selected loci as invalid for the other type
            if self._restr_sites:
                if is_qtl:
                    self.invalid_snp[chr_idx].extend(selected)
                else:
                    self.invalid_qtl[chr_idx].extend(selected)
        
        return loci
    
    def _samp_add_eff(self, qtl_loci: List[List[int]], n_traits: int, 
                     cor_a: np.ndarray, gamma: List[bool], shape: List[float]) -> np.ndarray:
        """
        Sample additive effects for QTL
        
        Parameters:
        -----------
        qtl_loci : list of list of int
            QTL loci for each chromosome
        n_traits : int
            Number of traits
        cor_a : numpy.ndarray
            Correlation matrix
        gamma : list of bool
            Whether to use gamma distribution for each trait
        shape : list of float
            Shape parameters for gamma distribution
        
        Returns:
        --------
        numpy.ndarray
            Additive effects matrix
        """
        n_qtl = sum(len(chr_loci) for chr_loci in qtl_loci)
        
        if n_traits == 1:
            if gamma[0]:
                # Gamma distribution
                effects = np.random.gamma(shape[0], 1, n_qtl)
            else:
                # Normal distribution
                effects = np.random.normal(0, 1, n_qtl)
            return effects.reshape(-1, 1)
        
        # Multiple traits with correlation
        # Generate correlated effects
        chol = np.linalg.cholesky(cor_a)
        effects = np.random.normal(0, 1, (n_qtl, n_traits))
        effects = effects @ chol.T
        
        # Apply gamma transformation if needed
        for i in range(n_traits):
            if gamma[i]:
                # Transform to gamma distribution
                effects[:, i] = np.random.gamma(shape[i], 1, n_qtl)
        
        return effects
    
    def _calc_gen_param(self, trait: TraitA, pop: MapPop) -> Dict[str, np.ndarray]:
        """
        Calculate genetic parameters for a trait
        
        Parameters:
        -----------
        trait : TraitA
            The trait to calculate parameters for
        pop : MapPop
            Population to calculate parameters on
        
        Returns:
        --------
        dict
            Dictionary containing breeding values and genetic values
        """
        # Calculate actual genetic values from QTL effects
        n_ind = pop.n_ind
        gv = np.zeros(n_ind)
        
        # Calculate genetic values for each individual
        for ind_idx in range(n_ind):
            total_gv = trait.intercept
            
            # Sum QTL effects across all QTL positions
            for qtl_idx, qtl_pos in enumerate(trait.qtl_loci):
                # For single chromosome, QTL positions are already local
                if pop.n_chr == 1:
                    local_pos = qtl_pos
                    chr_idx = 0
                else:
                    # Find which chromosome this QTL is on
                    chr_idx = 0
                    cumulative_loci = 0
                    for c in range(pop.n_chr):
                        if qtl_pos < cumulative_loci + pop.n_loci[c]:
                            chr_idx = c
                            break
                        cumulative_loci += pop.n_loci[c]
                    
                    # Get local position within chromosome
                    local_pos = qtl_pos - cumulative_loci
                
                if local_pos < pop.n_loci[chr_idx]:
                    # Get genotype at QTL position (unpack from binary format)
                    bin_idx = local_pos // 8
                    bit_idx = local_pos % 8
                    
                    # Extract genotype from packed format
                    dosage = 0
                    for ploidy_idx in range(pop.ploidy):
                        if pop.geno[chr_idx][bin_idx, ploidy_idx, ind_idx] & (1 << bit_idx):
                            dosage += 1
                    
                    # Calculate additive effect
                    # Scale dosage to [-1, 1] range for additive effects
                    scaled_dosage = (dosage - pop.ploidy/2) * (2/pop.ploidy)
                    
                    # Add QTL effect
                    if qtl_idx < len(trait.add_eff):
                        total_gv += trait.add_eff[qtl_idx] * scaled_dosage
            
            gv[ind_idx] = total_gv
        
        # For breeding values, use genetic values (simplified)
        bv = gv - trait.intercept
        
        return {'bv': bv, 'gv': gv}
    
    def _pop_var(self, x: np.ndarray) -> np.ndarray:
        """
        Calculate population variance (divide by n instead of n-1)
        
        Parameters:
        -----------
        x : numpy.ndarray
            Data array
        
        Returns:
        --------
        numpy.ndarray
            Population variance
        """
        return np.var(x, axis=0)
    
    def _add_trait(self, trait: TraitA, var_a: float, var_g: float):
        """
        Add trait to the simulation
        
        Parameters:
        -----------
        trait : TraitA
            Trait to add
        var_a : float
            Additive genetic variance
        var_g : float
            Total genetic variance
        """
        self._traits.append(trait)
        self._var_a.append(var_a)
        self._var_g.append(var_g)
        self._var_e.append(1.0)  # Default error variance


# Add convenience methods to SimParam class
def addTraitA(self, nQtlPerChr, mean=0, var=1, corA=None, gamma=False, shape=1, force=False, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_trait_a(nQtlPerChr, mean, var, corA, gamma, shape, force, name)

def addTraitAG(self, nQtlPerChr, mean=0, var=1, varGxE=1e-6, varEnv=0, corA=None, corGxE=None, gamma=False, shape=1, force=False, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_trait_ag(nQtlPerChr, mean, var, varGxE, varEnv, corA, corGxE, gamma, shape, force, name)

def restrSegSites(self, minQtlPerChr=None, minSnpPerChr=None, excludeQtl=None, excludeSnp=None, overlap=False, minSnpFreq=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.restr_seg_sites(minQtlPerChr, minSnpPerChr, excludeQtl, excludeSnp, overlap, minSnpFreq)

def addSnpChip(self, nSnpPerChr, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_snp_chip(nSnpPerChr, name)

def setTrackPed(self, isTrackPed, force=False):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.set_track_ped(isTrackPed, force)

def setTrackRec(self, isTrackRec, force=False):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.set_track_rec(isTrackRec, force)

def resetPed(self, lastId=0):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.reset_ped(lastId)

def setSexes(self, sexes, force=False):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.set_sexes(sexes, force)

# Add methods to SimParam class
SimParam.addTraitA = addTraitA
SimParam.addTraitAG = addTraitAG
SimParam.restrSegSites = restrSegSites
SimParam.addSnpChip = addSnpChip
SimParam.setTrackPed = setTrackPed
SimParam.setTrackRec = setTrackRec
SimParam.resetPed = resetPed
SimParam.setSexes = setSexes


def setVarE(self, h2=None, H2=None, varE=None, corE=None):
    """Set error variance from h2, H2, or varE. From AlphaSimR SimParam setVarE."""
    if varE is not None:
        varE = np.atleast_1d(varE)
        self._var_e = list(varE)
    elif h2 is not None:
        h2 = np.atleast_1d(h2)
        var_g = np.array(self._var_g)
        self._var_e = list(var_g * (1 - h2) / h2)
    elif H2 is not None:
        H2 = np.atleast_1d(H2)
        var_g = np.array(self._var_g)
        self._var_e = list(var_g * (1 - H2) / H2)
    else:
        raise ValueError("Must provide h2, H2, or varE")
    return self


SimParam.setVarE = setVarE


# Pop class
@dataclass
class Pop:
    """Population class extending MapPop with sex, genetic values, phenotypes, and pedigrees"""
    # Inherited from MapPop
    n_ind: int
    n_chr: int
    ploidy: int
    n_loci: List[int]
    geno: List[np.ndarray]
    gen_map: List[np.ndarray]
    centromere: List[float]
    inbred: bool
    
    # Additional Pop-specific attributes
    id: List[str]
    iid: List[int]
    mother: List[str]
    father: List[str]
    sex: List[str]
    n_traits: int
    gv: np.ndarray
    pheno: np.ndarray
    ebv: np.ndarray
    gxe: List[Any]
    fix_eff: List[int]
    misc: Dict[str, Any]
    misc_pop: Dict[str, Any]
    
    def __post_init__(self):
        """Validate the Pop object"""
        # Validate inherited MapPop attributes
        if self.n_chr != len(self.geno):
            raise ValueError("nChr != length(geno)")
        if self.n_chr != len(self.n_loci):
            raise ValueError("nChr != length(nLoci)")
        if self.n_chr != len(self.gen_map):
            raise ValueError("nChr != length(genMap)")
        if self.n_chr != len(self.centromere):
            raise ValueError("nChr != length(centromere)")
        
        # Validate Pop-specific attributes
        if self.n_ind != len(self.sex):
            raise ValueError("nInd != length(sex)")
        if self.n_ind != len(self.id):
            raise ValueError("nInd != length(id)")
        if self.n_ind != len(self.iid):
            raise ValueError("nInd != length(iid)")
        if self.n_ind != len(self.mother):
            raise ValueError("nInd != length(mother)")
        if self.n_ind != len(self.father):
            raise ValueError("nInd != length(father)")
        if self.n_ind != self.gv.shape[0]:
            raise ValueError("nInd != nrow(gv)")
        if self.n_ind != self.pheno.shape[0]:
            raise ValueError("nInd != nrow(pheno)")
        if self.n_ind != self.ebv.shape[0]:
            raise ValueError("nInd != nrow(ebv)")
        if self.n_traits != self.gv.shape[1]:
            raise ValueError("nTraits != ncol(gv)")
        if self.n_traits != self.pheno.shape[1]:
            raise ValueError("nTraits != ncol(pheno)")
        
        # Check for spaces in IDs
        for i, individual_id in enumerate(self.id):
            if ' ' in individual_id:
                raise ValueError(f"id[{i}] cannot contain spaces")
        for i, mother_id in enumerate(self.mother):
            if ' ' in mother_id:
                raise ValueError(f"mother[{i}] cannot contain spaces")
        for i, father_id in enumerate(self.father):
            if ' ' in father_id:
                raise ValueError(f"father[{i}] cannot contain spaces")
    
    def __getitem__(self, i):
        """Subset Pop by individual index or id. From AlphaSimR Pop [ method."""
        if isinstance(i, (slice, list, np.ndarray)):
            idx = np.arange(self.n_ind)[i]
        elif isinstance(i, str):
            idx = np.where(np.array(self.id) == i)[0]
            if len(idx) == 0:
                raise ValueError("Trying to select invalid individuals")
            idx = idx[0]
        else:
            idx = i
            if isinstance(idx, np.ndarray):
                idx = idx.tolist()
            if isinstance(idx, (list, np.ndarray)):
                idx = np.asarray(idx)
                if np.any(np.abs(idx) > self.n_ind):
                    raise ValueError("Trying to select invalid individuals")
                idx = (np.arange(self.n_ind)[idx] if np.any(idx < 0) else idx)
            else:
                if abs(idx) > self.n_ind:
                    raise ValueError("Trying to select invalid individuals")
                idx = np.arange(self.n_ind)[idx] if idx < 0 else idx
        idx = np.atleast_1d(idx)
        if idx.ndim == 0:
            idx = idx.reshape(1)
        subset_misc = {}
        for k, v in self.misc.items():
            if isinstance(v, np.ndarray) and v.ndim > 0:
                subset_misc[k] = v[idx] if v.shape[0] == self.n_ind else v
            else:
                subset_misc[k] = [v[i] for i in idx] if hasattr(v, '__getitem__') else v
        return Pop(
            n_ind=len(idx), n_chr=self.n_chr, ploidy=self.ploidy, n_loci=self.n_loci,
            geno=[g[:, :, idx] for g in self.geno], gen_map=self.gen_map,
            centromere=self.centromere, inbred=self.inbred,
            id=[self.id[j] for j in idx], iid=[self.iid[j] for j in idx],
            mother=[self.mother[j] for j in idx], father=[self.father[j] for j in idx],
            sex=[self.sex[j] for j in idx], n_traits=self.n_traits,
            gv=self.gv[idx], pheno=self.pheno[idx], ebv=self.ebv[idx],
            gxe=self.gxe, fix_eff=[self.fix_eff[j] for j in idx],
            misc=subset_misc, misc_pop={}
        )


def merge_pops(pop_list) -> Pop:
    """Merge list of Pop objects. From AlphaSimR mergePops.R"""
    if hasattr(pop_list, 'pops'):
        pop_list = pop_list.pops
    pop_list = [p for p in pop_list if p is not None and isinstance(p, Pop)]
    if not pop_list:
        raise ValueError("popList must contain Pop objects")
    n_chr = pop_list[0].n_chr
    ploidy = pop_list[0].ploidy
    n_loci = pop_list[0].n_loci
    for p in pop_list[1:]:
        if p.n_chr != n_chr or p.ploidy != ploidy or not all(np.array(p.n_loci) == np.array(n_loci)):
            raise ValueError("All populations must have compatible structure")
    misc = {}
    if all(len(p.misc) == len(pop_list[0].misc) for p in pop_list) and len(pop_list[0].misc) > 0:
        if all(set(p.misc.keys()) == set(pop_list[0].misc.keys()) for p in pop_list):
            for k in pop_list[0].misc.keys():
                vals = [p.misc[k] for p in pop_list]
                if isinstance(vals[0], np.ndarray) and vals[0].ndim > 0:
                    misc[k] = np.vstack(vals) if vals[0].ndim > 1 else np.concatenate(vals)
                else:
                    misc[k] = np.concatenate([np.atleast_1d(v) for v in vals])
    n_ind = sum(p.n_ind for p in pop_list)
    geno = []
    for chr_idx in range(n_chr):
        geno.append(np.concatenate([p.geno[chr_idx] for p in pop_list], axis=2))
    return Pop(
        n_ind=n_ind, n_chr=n_chr, ploidy=ploidy, n_loci=n_loci,
        geno=geno, gen_map=pop_list[0].gen_map, centromere=pop_list[0].centromere,
        inbred=pop_list[0].inbred,
        id=[x for p in pop_list for x in p.id],
        iid=[x for p in pop_list for x in p.iid],
        mother=[x for p in pop_list for x in p.mother],
        father=[x for p in pop_list for x in p.father],
        sex=[x for p in pop_list for x in p.sex],
        n_traits=pop_list[0].n_traits,
        gv=np.vstack([p.gv for p in pop_list]),
        pheno=np.vstack([p.pheno for p in pop_list]),
        ebv=np.vstack([p.ebv for p in pop_list]) if pop_list[0].ebv.shape[1] > 0 else np.empty((n_ind, 0)),
        gxe=pop_list[0].gxe,
        fix_eff=[x for p in pop_list for x in p.fix_eff],
        misc=misc, misc_pop={}
    )


def new_pop(raw_pop: MapPop, sim_param: Optional[SimParam] = None, **kwargs) -> Pop:
    """
    Create new population from MapPop
    
    Parameters:
    -----------
    raw_pop : MapPop
        The raw population (MapPop) to convert
    sim_param : SimParam, optional
        Simulation parameters. If None, uses global SP
    **kwargs
        Additional arguments passed to finalizePop function
    
    Returns:
    --------
    Pop
        New population object
    """
    if sim_param is None:
        # In a full implementation, this would look for global SP
        raise ValueError("simParam must be provided")
    
    # Validate genetic map compatibility
    for i in range(raw_pop.n_chr):
        if len(sim_param._female_map[i]) != raw_pop.n_loci[i]:
            raise ValueError(f"Genetic map length mismatch for chromosome {i+1}")
    
    # Generate IDs
    last_id = sim_param._last_id
    iid = list(range(last_id + 1, last_id + raw_pop.n_ind + 1))
    last_id = max(iid)
    
    # Generate individual IDs
    id_list = [str(i) for i in iid]
    
    # Generate parent IDs (all founders have parents "0")
    mother = ["0"] * raw_pop.n_ind
    father = ["0"] * raw_pop.n_ind
    
    # Generate internal parent IDs
    i_mother = [0] * raw_pop.n_ind
    i_father = [0] * raw_pop.n_ind
    
    # Determine if individuals are DH/inbred
    is_dh = raw_pop.inbred
    
    # Generate sex
    if sim_param._sexes == "no":
        sex = ["H"] * raw_pop.n_ind
    elif sim_param._sexes == "yes_rand":
        sex = np.random.choice(["M", "F"], raw_pop.n_ind, replace=True).tolist()
    elif sim_param._sexes == "yes_sys":
        sex = (["M", "F"] * ((raw_pop.n_ind // 2) + 1))[:raw_pop.n_ind]
    else:
        raise ValueError(f"No rules for sex type {sim_param._sexes}")
    
    # Initialize genetic values and phenotypes
    n_traits = sim_param.n_traits
    gv = np.full((raw_pop.n_ind, n_traits), np.nan)
    pheno = np.full((raw_pop.n_ind, n_traits), np.nan)
    ebv = np.empty((raw_pop.n_ind, 0))
    gxe = [None] * n_traits
    
    # Calculate genetic values if traits exist
    if n_traits >= 1:
        gv = _get_gv_index(raw_pop, sim_param)
        
        # Set phenotypes (simplified - just genetic values for now)
        pheno = gv.copy()
    
    # Create Pop object
    pop = Pop(
        # Inherited from MapPop
        n_ind=raw_pop.n_ind,
        n_chr=raw_pop.n_chr,
        ploidy=raw_pop.ploidy,
        n_loci=raw_pop.n_loci,
        geno=raw_pop.geno,
        gen_map=raw_pop.gen_map,
        centromere=raw_pop.centromere,
        inbred=raw_pop.inbred,
        
        # Pop-specific attributes
        id=id_list,
        iid=iid,
        mother=mother,
        father=father,
        sex=sex,
        n_traits=n_traits,
        gv=gv,
        pheno=pheno,
        ebv=ebv,
        gxe=gxe,
        fix_eff=[1] * raw_pop.n_ind,
        misc={},
        misc_pop={}
    )
    
    # Apply finalizePop function
    pop = sim_param.finalize_pop(pop)
    
    # Update tracking if enabled
    if sim_param._is_track_ped:
        if sim_param._is_track_rec:
            # Add to recombination history (simplified)
            sim_param._rec_hist.extend([None] * raw_pop.n_ind)
        else:
            # Add to pedigree (simplified)
            pass
    
    # Update last ID
    sim_param._last_id = last_id
    
    return pop


def _get_gv_index(raw_pop: MapPop, sim_param: SimParam) -> np.ndarray:
    """
    Calculate genetic values for a population based on QTL effects and genotypes
    
    Parameters:
    -----------
    raw_pop : MapPop
        The raw population
    sim_param : SimParam
        Simulation parameters
    
    Returns:
    --------
    numpy.ndarray
        Genetic values matrix (n_ind x n_traits)
    """
    n_ind = raw_pop.n_ind
    n_traits = sim_param.n_traits
    
    if n_traits == 0:
        return np.empty((n_ind, 0))
    
    gv = np.zeros((n_ind, n_traits))
    
    # Calculate genetic values based on QTL effects and genotypes
    for trait_idx in range(n_traits):
        trait = sim_param._traits[trait_idx]
        
        # Get QTL loci for this trait (flat list of positions)
        qtl_loci = trait.qtl_loci
        
        # Calculate genetic values for each individual
        for ind_idx in range(n_ind):
            total_gv = trait.intercept
            
            # Sum QTL effects across all QTL positions
            for qtl_idx, qtl_pos in enumerate(qtl_loci):
                # For single chromosome, QTL positions are already local
                if raw_pop.n_chr == 1:
                    local_pos = qtl_pos
                    chr_idx = 0
                else:
                    # Find which chromosome this QTL is on
                    chr_idx = 0
                    cumulative_loci = 0
                    for c in range(raw_pop.n_chr):
                        if qtl_pos < cumulative_loci + raw_pop.n_loci[c]:
                            chr_idx = c
                            break
                        cumulative_loci += raw_pop.n_loci[c]
                    
                    # Get local position within chromosome
                    local_pos = qtl_pos - cumulative_loci
                
                if local_pos < raw_pop.n_loci[chr_idx]:
                    # Get genotype at QTL position (unpack from binary format)
                    bin_idx = local_pos // 8
                    bit_idx = local_pos % 8
                    
                    # Extract genotype from packed format
                    dosage = 0
                    for ploidy_idx in range(raw_pop.ploidy):
                        if raw_pop.geno[chr_idx][bin_idx, ploidy_idx, ind_idx] & (1 << bit_idx):
                            dosage += 1
                    
                    # Calculate additive effect
                    # Scale dosage to [-1, 1] range for additive effects
                    scaled_dosage = (dosage - raw_pop.ploidy/2) * (2/raw_pop.ploidy)
                    
                    # Add QTL effect
                    if qtl_idx < len(trait.add_eff):
                        total_gv += trait.add_eff[qtl_idx] * scaled_dosage
            
            gv[ind_idx, trait_idx] = total_gv
    
    return gv


# Genetic statistics functions
def mean_g(pop: Pop) -> np.ndarray:
    """
    Calculate mean genetic values for all traits
    
    Parameters:
    -----------
    pop : Pop
        Population object
    
    Returns:
    --------
    numpy.ndarray
        Mean genetic values for each trait
    """
    return np.mean(pop.gv, axis=0)


def var_g(pop: Pop) -> np.ndarray:
    """
    Calculate genetic variance for all traits
    
    Parameters:
    -----------
    pop : Pop
        Population object
    
    Returns:
    --------
    numpy.ndarray
        Genetic variance for each trait
    """
    return np.var(pop.gv, axis=0)


def mean_p(pop: Pop) -> np.ndarray:
    """
    Calculate mean phenotypic values for all traits
    
    Parameters:
    -----------
    pop : Pop
        Population object
    
    Returns:
    --------
    numpy.ndarray
        Mean phenotypic values for each trait
    """
    return np.mean(pop.pheno, axis=0)


# Phenotype functions
def set_pheno(pop: Pop, var_e: Optional[Union[float, List[float], np.ndarray]] = None,
              reps: int = 1, sim_param: Optional[SimParam] = None) -> Pop:
    """
    Set phenotypes for all traits by adding random error
    
    Parameters:
    -----------
    pop : Pop
        Population object
    var_e : float, list of float, or numpy.ndarray, optional
        Error variance for each trait
    reps : int, default=1
        Number of replications for phenotype
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop
        Population with updated phenotypes
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    n_traits = pop.n_traits
    if n_traits == 0:
        return pop
    
    # Set error variance
    if var_e is None:
        var_e = sim_param._var_e
    elif isinstance(var_e, (int, float)):
        var_e = [var_e] * n_traits
    elif isinstance(var_e, np.ndarray):
        if var_e.ndim == 1:
            var_e = var_e.tolist()
        else:
            # Covariance matrix - extract diagonal
            var_e = np.diag(var_e).tolist()
    
    # Generate phenotypes
    n_ind = pop.n_ind
    pheno = np.zeros((n_ind, n_traits))
    
    for trait_idx in range(n_traits):
        # Genetic value + error
        error_var = var_e[trait_idx] / reps  # Scale by replications
        error = np.random.normal(0, np.sqrt(error_var), n_ind)
        pheno[:, trait_idx] = pop.gv[:, trait_idx] + error
    
    # Update population
    pop.pheno = pheno
    return pop


def _cut_r_style(x: np.ndarray, breaks: np.ndarray, include_lowest: bool, right: bool) -> np.ndarray:
    """Replicate R's cut() behavior. Returns 1-indexed categories, NaN for out of range."""
    result = np.full_like(x, np.nan, dtype=float)
    breaks = np.asarray(breaks)
    n_breaks = len(breaks)
    if n_breaks < 2:
        return result
    # right=FALSE: [a,b) intervals. include.lowest=TRUE: last interval includes upper bound.
    for i in range(n_breaks - 1):
        lo, hi = breaks[i], breaks[i + 1]
        if right:
            # (a,b] - exclude left, include right
            if i == 0 and include_lowest:
                mask = (x >= lo) & (x <= hi)
            else:
                mask = (x > lo) & (x <= hi)
        else:
            # [a,b) - include left, exclude right
            if i == n_breaks - 2 and include_lowest:
                mask = (x >= lo) & (x <= hi)
            else:
                mask = (x >= lo) & (x < hi)
        result[mask] = i + 1
    return result


def as_categorical(x: np.ndarray, p: Optional[Union[float, np.ndarray, List]] = None,
                   mean: Union[float, np.ndarray] = 0, var: Union[float, np.ndarray] = 1,
                   threshold: Union[np.ndarray, List] = None,
                   include_lowest: bool = True, right: bool = False) -> np.ndarray:
    """
    Convert continuous (Gaussian) trait to categorical via ordered probit.
    From AlphaSimR phenotypes.R asCategorical.
    """
    import warnings
    from scipy import stats
    
    x = np.asarray(x)
    if x.ndim == 1:
        x = x.reshape(-1, 1)
    n_traits = x.shape[1]
    
    if p is not None:
        if np.isscalar(p):
            if n_traits > 1:
                raise ValueError("When x contains more than one column, you must supply a list of probabilities!")
            if p == 0.5 and n_traits == 1:
                p = [np.array([0.5, 0.5])]
            else:
                p = [p]
        elif isinstance(p, (list, tuple)) and len(p) == n_traits:
            p = list(p)
        else:
            p_arr = np.atleast_1d(p)
            # When n_traits==1 and p is category probs (e.g. [0.5, 0.5]), wrap as [p]
            if n_traits == 1 and len(p_arr) > 0 and np.isscalar(p_arr.flat[0]):
                p = [p_arr]
            else:
                p = list(p) if hasattr(p, '__iter__') and not isinstance(p, np.ndarray) else [p]
        if len(p) != n_traits:
            raise ValueError("You must supply probabilities for all traits in x!")
        mean_arr = np.atleast_1d(mean)
        var_arr = np.atleast_1d(var)
        if any(pp is not None for pp in p):
            if len(mean_arr) != n_traits:
                raise ValueError("You must supply means for all traits in x!")
            if len(var_arr) != n_traits:
                raise ValueError("You must supply variances for all traits in x!")
        threshold = []
        for trt in range(n_traits):
            p_trt = p[trt]
            if p_trt is None:
                threshold.append(None)
                continue
            p_trt = np.atleast_1d(p_trt)
            p_sum = np.sum(p_trt)
            if not np.isclose(p_sum, 1):
                warnings.warn("Probabilities do not sum to 1 - creating one more category!")
                p_trt = np.concatenate([p_trt, [1 - p_sum]])
            tmp = stats.norm.ppf(np.cumsum(p_trt), loc=mean_arr[trt], scale=np.sqrt(var_arr[trt]))
            if not np.any(np.isinf(tmp) & (tmp < 0)):
                tmp = np.concatenate([[-np.inf], tmp])
            if not np.any(np.isinf(tmp) & (tmp > 0)):
                tmp = np.concatenate([tmp, [np.inf]])
            threshold.append(tmp)
    
    if threshold is None:
        if n_traits > 1:
            raise ValueError("When x contains more than one column, you must supply a list of thresholds!")
        threshold = [np.array([-np.inf, 0, np.inf])]
    
    if not isinstance(threshold, (list, tuple)):
        threshold = [np.asarray(threshold)]
    else:
        threshold = [np.asarray(t) if t is not None else None for t in threshold]
    
    if len(threshold) != n_traits:
        raise ValueError("You must supply thresholds for all traits in x!")
    
    result = x.copy().astype(float)
    for trt in range(n_traits):
        if threshold[trt] is not None:
            result[:, trt] = _cut_r_style(
                x[:, trt], threshold[trt], include_lowest, right
            )
    return result


# Selection functions
def select_ind(pop: Pop, n_ind: int, trait: Union[int, callable] = 1,
               use: Union[str, callable] = "pheno", sex: str = "B",
               select_top: bool = True, return_pop: bool = True,
               candidates: Optional[List[int]] = None,
               parents: Optional[List[int]] = None,
               sim_param: Optional[SimParam] = None, **kwargs) -> Union[Pop, List[int]]:
    """
    Select a subset of individuals from a population
    
    Parameters:
    -----------
    pop : Pop
        Population object
    n_ind : int
        Number of individuals to select
    trait : int or callable, default=1
        Trait for selection (trait index or function)
    use : str or callable, default="pheno"
        Selection criterion ("gv", "pheno", "ebv", "rand")
    sex : str, default="B"
        Which sex to select ("B", "F", "M")
    select_top : bool, default=True
        Select highest values if True
    return_pop : bool, default=True
        Return Pop object if True, indices if False
    candidates : list of int, optional
        Eligible selection candidates
    sim_param : SimParam, optional
        Simulation parameters
    **kwargs
        Additional arguments for trait/use functions
    
    Returns:
    --------
    Pop or list of int
        Selected population or indices
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")

    if candidates is None and parents is not None:
        candidates = parents

    if n_ind <= 0:
        if return_pop:
            # Return empty population
            return Pop(
                n_ind=0, n_chr=pop.n_chr, ploidy=pop.ploidy, n_loci=pop.n_loci,
                geno=[], gen_map=pop.gen_map, centromere=pop.centromere,
                inbred=pop.inbred, id=[], iid=[], mother=[], father=[],
                sex=[], n_traits=pop.n_traits, gv=np.empty((0, pop.n_traits)),
                pheno=np.empty((0, pop.n_traits)), ebv=np.empty((0, 0)),
                gxe=[], fix_eff=[], misc={}, misc_pop={}
            )
        else:
            return []
    
    # Get eligible individuals based on sex
    eligible = []
    for i in range(pop.n_ind):
        if sex == "B" or pop.sex[i] == sex:
            eligible.append(i)
    
    # Filter by candidates if provided
    if candidates is not None:
        eligible = [i for i in eligible if i in candidates]
    
    if len(eligible) < n_ind:
        n_ind = len(eligible)
        if n_ind == 0:
            if return_pop:
                return Pop(
                    n_ind=0, n_chr=pop.n_chr, ploidy=pop.ploidy, n_loci=pop.n_loci,
                    geno=[], gen_map=pop.gen_map, centromere=pop.centromere,
                    inbred=pop.inbred, id=[], iid=[], mother=[], father=[],
                    sex=[], n_traits=pop.n_traits, gv=np.empty((0, pop.n_traits)),
                    pheno=np.empty((0, pop.n_traits)), ebv=np.empty((0, 0)),
                    gxe=[], fix_eff=[], misc={}, misc_pop={}
                )
            else:
                return []
    
    # Get selection criterion values
    if isinstance(use, str):
        if use == "pheno":
            values = pop.pheno[:, trait-1] if isinstance(trait, int) else pop.pheno
        elif use == "gv":
            values = pop.gv[:, trait-1] if isinstance(trait, int) else pop.gv
        elif use == "ebv":
            values = pop.ebv[:, trait-1] if isinstance(trait, int) else pop.ebv
        elif use == "rand":
            values = np.random.random(pop.n_ind)
        else:
            raise ValueError(f"Unknown use criterion: {use}")
    else:
        # Use function
        values = use(pop, trait, **kwargs)
    
    # Apply trait function if needed
    if callable(trait):
        values = trait(values, **kwargs)
    
    # Ensure values is 1D
    if values.ndim > 1:
        values = values.flatten()
    
    # Select individuals (match R order(): decreasing=selectTop gives highest first)
    eligible_values = values[eligible]
    sorted_indices = np.argsort(eligible_values)
    if select_top:
        selected_indices = sorted_indices[-n_ind:][::-1]
    else:
        selected_indices = sorted_indices[:n_ind]
    
    selected_eligible = [eligible[i] for i in selected_indices]
    
    if return_pop:
        # Create new population with selected individuals
        return Pop(
            n_ind=n_ind, n_chr=pop.n_chr, ploidy=pop.ploidy, n_loci=pop.n_loci,
            geno=[pop.geno[chr_idx][:, :, selected_eligible] for chr_idx in range(pop.n_chr)],
            gen_map=pop.gen_map, centromere=pop.centromere, inbred=pop.inbred,
            id=[pop.id[i] for i in selected_eligible],
            iid=[pop.iid[i] for i in selected_eligible],
            mother=[pop.mother[i] for i in selected_eligible],
            father=[pop.father[i] for i in selected_eligible],
            sex=[pop.sex[i] for i in selected_eligible],
            n_traits=pop.n_traits, gv=pop.gv[selected_eligible, :],
            pheno=pop.pheno[selected_eligible, :], ebv=pop.ebv[selected_eligible, :],
            gxe=pop.gxe, fix_eff=[pop.fix_eff[i] for i in selected_eligible],
            misc=pop.misc, misc_pop=pop.misc_pop
        )
    else:
        return selected_eligible


# Crossing functions
def rand_cross(pop: Pop, n_crosses: int, n_progeny: int = 1,
               balance: bool = True, parents: Optional[List[int]] = None,
               ignore_sexes: bool = False, sim_param: Optional[SimParam] = None) -> Pop:
    """
    Make random crosses between individuals in a population
    
    Parameters:
    -----------
    pop : Pop
        Population object
    n_crosses : int
        Number of crosses to make
    n_progeny : int, default=1
        Number of progeny per cross
    balance : bool, default=True
        Balance number of progeny per parent
    parents : list of int, optional
        Allowable parent indices
    ignore_sexes : bool, default=False
        Ignore sex when making crosses
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop
        New population of progeny
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    n_ind = pop.n_ind
    if n_ind < 2:
        raise ValueError("Need at least 2 individuals to make crosses")
    
    # Determine eligible parents
    if parents is None:
        eligible_parents = list(range(n_ind))
    else:
        eligible_parents = parents
    
    if len(eligible_parents) < 2:
        raise ValueError("Need at least 2 eligible parents")
    
    # Separate by sex if not ignoring sexes
    if not ignore_sexes and sim_param._sexes != "no":
        females = [i for i in eligible_parents if pop.sex[i] == "F"]
        males = [i for i in eligible_parents if pop.sex[i] == "M"]
        
        if len(females) == 0 or len(males) == 0:
            raise ValueError("Need both male and female parents when using sexes")
        
        # Make crosses
        cross_plan = []
        for _ in range(n_crosses):
            female_idx = np.random.choice(females)
            male_idx = np.random.choice(males)
            cross_plan.append((female_idx, male_idx))
    else:
        # Random crosses ignoring sex
        cross_plan = []
        for _ in range(n_crosses):
            parents_idx = np.random.choice(eligible_parents, size=2, replace=False)
            cross_plan.append((parents_idx[0], parents_idx[1]))
    
    # Generate progeny
    total_progeny = n_crosses * n_progeny
    progeny_geno = []
    
    for chr_idx in range(pop.n_chr):
        # Get the packed genotype dimensions
        n_packed = pop.geno[chr_idx].shape[0]
        chr_geno = np.zeros((n_packed, pop.ploidy, total_progeny), dtype=np.uint8)
        
        for cross_idx, (female_idx, male_idx) in enumerate(cross_plan):
            for prog_idx in range(n_progeny):
                prog_idx_global = cross_idx * n_progeny + prog_idx
                
                # Simple crossing simulation - randomly inherit from each parent
                for packed_idx in range(n_packed):
                    for haplo in range(pop.ploidy):
                        if np.random.random() < 0.5:
                            # Inherit from female
                            chr_geno[packed_idx, haplo, prog_idx_global] = pop.geno[chr_idx][packed_idx, haplo, female_idx]
                        else:
                            # Inherit from male
                            chr_geno[packed_idx, haplo, prog_idx_global] = pop.geno[chr_idx][packed_idx, haplo, male_idx]
        
        progeny_geno.append(chr_geno)
    
    # Create progeny population
    progeny_pop = Pop(
        n_ind=total_progeny, n_chr=pop.n_chr, ploidy=pop.ploidy, n_loci=pop.n_loci,
        geno=progeny_geno, gen_map=pop.gen_map, centromere=pop.centromere,
        inbred=pop.inbred,
        id=[f"P{i+1}" for i in range(total_progeny)],
        iid=list(range(sim_param._last_id + 1, sim_param._last_id + total_progeny + 1)),
        mother=[pop.id[female_idx] for female_idx, _ in cross_plan for _ in range(n_progeny)],
        father=[pop.id[male_idx] for _, male_idx in cross_plan for _ in range(n_progeny)],
        sex=["H"] * total_progeny,  # Default to hermaphrodite
        n_traits=pop.n_traits,
        gv=np.zeros((total_progeny, pop.n_traits)),
        pheno=np.zeros((total_progeny, pop.n_traits)),
        ebv=np.empty((total_progeny, 0)),
        gxe=[None] * pop.n_traits,
        fix_eff=[1] * total_progeny,
        misc={},
        misc_pop={}
    )
    
    # Calculate genetic values for progeny
    if pop.n_traits > 0:
        progeny_pop.gv = _get_gv_index(progeny_pop, sim_param)
        progeny_pop.pheno = progeny_pop.gv.copy()
    
    # Update last ID
    sim_param._last_id += total_progeny
    
    return progeny_pop


def _do_crosses(females: Pop, males: Pop, cross_plan: np.ndarray, n_progeny: int,
                 sim_param: SimParam) -> Pop:
    """Execute crosses from cross plan. Returns progeny Pop."""
    total_progeny = len(cross_plan) * n_progeny
    progeny_geno = []
    for chr_idx in range(females.n_chr):
        n_packed = females.geno[chr_idx].shape[0]
        chr_geno = np.zeros((n_packed, females.ploidy, total_progeny), dtype=np.uint8)
        for cross_idx, (f_idx, m_idx) in enumerate(cross_plan):
            for prog_idx in range(n_progeny):
                prog_global = cross_idx * n_progeny + prog_idx
                for h in range(females.ploidy):
                    if h % 2 == 0:
                        hp = np.random.randint(0, females.ploidy)
                        chr_geno[:, h, prog_global] = females.geno[chr_idx][:, hp, f_idx]
                    else:
                        hp = np.random.randint(0, males.ploidy)
                        chr_geno[:, h, prog_global] = males.geno[chr_idx][:, hp, m_idx]
        progeny_geno.append(chr_geno)
    progeny_pop = Pop(
        n_ind=total_progeny, n_chr=females.n_chr, ploidy=females.ploidy, n_loci=females.n_loci,
        geno=progeny_geno, gen_map=females.gen_map, centromere=females.centromere,
        inbred=females.inbred,
        id=[f"P{i+1}" for i in range(total_progeny)],
        iid=list(range(sim_param._last_id + 1, sim_param._last_id + total_progeny + 1)),
        mother=[females.id[f] for f, _ in cross_plan for _ in range(n_progeny)],
        father=[males.id[m] for _, m in cross_plan for _ in range(n_progeny)],
        sex=["H"] * total_progeny, n_traits=females.n_traits,
        gv=np.zeros((total_progeny, females.n_traits)),
        pheno=np.zeros((total_progeny, females.n_traits)),
        ebv=np.empty((total_progeny, 0)), gxe=[None] * females.n_traits,
        fix_eff=[1] * total_progeny, misc={}, misc_pop={}
    )
    if females.n_traits > 0:
        progeny_pop.gv = _get_gv_index(progeny_pop, sim_param)
        progeny_pop.pheno = progeny_pop.gv.copy()
    sim_param._last_id += total_progeny
    return progeny_pop


def make_cross(pop: Pop, cross_plan: np.ndarray, n_progeny: int = 1,
               sim_param: Optional[SimParam] = None) -> Pop:
    """Make designed crosses. From AlphaSimR crossing.R makeCross."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    cross_plan = np.asarray(cross_plan)
    if cross_plan.dtype.kind in ('U', 'S', 'O'):
        def _lookup(s, id_list, n_ind):
            s = str(s)
            try:
                idx = int(s)
                if 0 <= idx < n_ind:
                    return idx
            except (ValueError, TypeError):
                pass
            try:
                return id_list.index(s)
            except ValueError:
                raise ValueError(f"'{s}' is not in list")
        f_idx = np.array([_lookup(r[0], pop.id, pop.n_ind) for r in cross_plan])
        m_idx = np.array([_lookup(r[1], pop.id, pop.n_ind) for r in cross_plan])
        cross_plan = np.column_stack([f_idx, m_idx])
    else:
        cross_plan = cross_plan.astype(int)
        if np.min(cross_plan) >= 1:
            cross_plan = cross_plan - 1
        if np.any(cross_plan < 0) or np.any(cross_plan >= pop.n_ind):
            raise ValueError("Invalid crossPlan")
    if n_progeny > 1:
        cross_plan = np.repeat(cross_plan, n_progeny, axis=0)
    return _do_crosses(pop, pop, cross_plan, 1, sim_param)


def make_cross2(females: Pop, males: Pop, cross_plan: np.ndarray, n_progeny: int = 1,
                sim_param: Optional[SimParam] = None) -> Pop:
    """Make crosses between two populations. From AlphaSimR crossing.R makeCross2."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    cross_plan = np.asarray(cross_plan)
    if cross_plan.dtype.kind in ('U', 'S', 'O'):
        def _lookup(s, id_list, n_ind):
            s = str(s)
            try:
                idx = int(s)
                if 0 <= idx < n_ind:
                    return idx
            except (ValueError, TypeError):
                pass
            try:
                return id_list.index(s)
            except ValueError:
                raise ValueError(f"'{s}' is not in list")
        f_idx = np.array([_lookup(r[0], females.id, females.n_ind) for r in cross_plan])
        m_idx = np.array([_lookup(r[1], males.id, males.n_ind) for r in cross_plan])
        cross_plan = np.column_stack([f_idx, m_idx])
    else:
        cross_plan = cross_plan.astype(int)
        if np.min(cross_plan) >= 1:
            cross_plan = cross_plan - 1
        if np.any(cross_plan[:, 0] < 0) or np.any(cross_plan[:, 0] >= females.n_ind):
            raise ValueError("Invalid crossPlan")
        if np.any(cross_plan[:, 1] < 0) or np.any(cross_plan[:, 1] >= males.n_ind):
            raise ValueError("Invalid crossPlan")
    if n_progeny > 1:
        cross_plan = np.repeat(cross_plan, n_progeny, axis=0)
    return _do_crosses(females, males, cross_plan, 1, sim_param)


def rand_cross2(females: Pop, males: Pop, n_crosses: int, n_progeny: int = 1,
                balance: bool = True, female_parents=None, male_parents=None,
                ignore_sexes: bool = False, sim_param: Optional[SimParam] = None) -> Pop:
    """Random crosses between two populations. From AlphaSimR crossing.R randCross2."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    female_parents = female_parents or list(range(females.n_ind))
    male_parents = male_parents or list(range(males.n_ind))
    if sim_param._sexes == "no" or ignore_sexes:
        f_eligible = female_parents
        m_eligible = male_parents
    else:
        f_eligible = [i for i in female_parents if females.sex[i] == "F"]
        m_eligible = [i for i in male_parents if males.sex[i] == "M"]
        if not f_eligible or not m_eligible:
            raise ValueError("Need both female and male parents")
    if balance:
        f_eligible = list(np.random.permutation(f_eligible))
        m_eligible = list(np.random.permutation(m_eligible))
        f_eligible = (f_eligible * (n_crosses // len(f_eligible) + 1))[:n_crosses]
        m_eligible = (m_eligible * (n_crosses // len(m_eligible) + 1))[:n_crosses]
        cross_plan = np.column_stack([f_eligible, m_eligible])
    else:
        cross_plan = np.column_stack([
            np.random.choice(f_eligible, n_crosses),
            np.random.choice(m_eligible, n_crosses)
        ])
    return _do_crosses(females, males, cross_plan, n_progeny, sim_param)


def self_(pop: Pop, n_progeny: int = 1, parents=None, keep_parents: bool = True,
          sim_param: Optional[SimParam] = None) -> Pop:
    """Self individuals. From AlphaSimR crossing.R self."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    parents = parents or list(range(pop.n_ind))
    parents = np.atleast_1d(parents).astype(int) - 1
    cross_plan = np.column_stack([np.repeat(parents, n_progeny),
                                  np.repeat(parents, n_progeny)])
    return _do_crosses(pop, pop, cross_plan, 1, sim_param)


def make_dh(pop: Pop, n_dh: int = 1, use_female: bool = True, keep_parents: bool = True,
            sim_param: Optional[SimParam] = None) -> Pop:
    """Create doubled haploids. From AlphaSimR crossing.R makeDH."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    if pop.ploidy != 2:
        raise ValueError("Only works with diploids")
    result_pops = []
    for ind_idx in range(pop.n_ind):
        for _ in range(n_dh):
            haplo_idx = 0 if np.random.random() < 0.5 else 1
            progeny_geno = []
            for chr_idx in range(pop.n_chr):
                g = pop.geno[chr_idx][:, :, ind_idx]
                dh_geno = np.stack([g[:, haplo_idx].copy(), g[:, haplo_idx].copy()], axis=1)
                dh_geno = dh_geno.reshape(dh_geno.shape[0], 2, 1)
                progeny_geno.append(dh_geno)
            result_pops.append(progeny_geno)
    total = pop.n_ind * n_dh
    merged_geno = [np.concatenate([p[c] for p in result_pops], axis=2) for c in range(pop.n_chr)]
    progeny_pop = Pop(
        n_ind=total, n_chr=pop.n_chr, ploidy=2, n_loci=pop.n_loci,
        geno=merged_geno, gen_map=pop.gen_map, centromere=pop.centromere,
        inbred=True,
        id=[f"P{i+1}" for i in range(total)],
        iid=list(range(sim_param._last_id + 1, sim_param._last_id + total + 1)),
        mother=[pop.mother[i] for i in range(pop.n_ind) for _ in range(n_dh)],
        father=[pop.father[i] for i in range(pop.n_ind) for _ in range(n_dh)],
        sex=["H"] * total, n_traits=pop.n_traits,
        gv=np.zeros((total, pop.n_traits)), pheno=np.zeros((total, pop.n_traits)),
        ebv=np.empty((total, 0)), gxe=[None] * pop.n_traits,
        fix_eff=[1] * total, misc={}, misc_pop={}
    )
    if pop.n_traits > 0:
        progeny_pop.gv = _get_gv_index(progeny_pop, sim_param)
        progeny_pop.pheno = progeny_pop.gv.copy()
    sim_param._last_id += total
    return progeny_pop


def select_cross(pop: Pop, n_ind=None, n_female=None, n_male=None, n_crosses: int = 1,
                n_progeny: int = 1, trait=1, use: str = "pheno", select_top: bool = True,
                balance: bool = True, sim_param: Optional[SimParam] = None, **kwargs) -> Pop:
    """Select and cross. From AlphaSimR crossing.R selectCross."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    if n_ind is not None:
        parents = select_ind(pop, n_ind, trait, use, "B", select_top, False, sim_param=sim_param, **kwargs)
        cross_plan = []
        for _ in range(n_crosses):
            p = np.random.choice(parents, 2, replace=False)
            cross_plan.append(p + 1)
        cross_plan = np.array(cross_plan)
        return make_cross(pop, cross_plan, n_progeny, sim_param)
    else:
        if sim_param._sexes == "no":
            raise ValueError("Must specify nInd when simParam sexes is 'no'")
        females = select_ind(pop, n_female, trait, use, "F", select_top, True, sim_param=sim_param, **kwargs)
        males = select_ind(pop, n_male, trait, use, "M", select_top, True, sim_param=sim_param, **kwargs)
        return rand_cross2(females, males, n_crosses, n_progeny, balance, sim_param=sim_param)


def select_op(pop: Pop, n_ind: int, n_seeds: int, prob_self: float = 0,
              pollen_control: bool = False, trait=1, use: str = "pheno",
              select_top: bool = True, candidates=None, sim_param: Optional[SimParam] = None, **kwargs) -> Pop:
    """Open pollination. From AlphaSimR selection.R selectOP."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    female = select_ind(
        pop,
        n_ind,
        trait,
        use,
        "B",
        select_top,
        False,
        candidates=candidates,
        sim_param=sim_param,
        **kwargs,
    )
    n_self = np.random.binomial(n_seeds, prob_self, n_ind)
    cross_plan = []
    for i, f in enumerate(female):
        males = [j for j in range(pop.n_ind) if j != f] if not pollen_control else [f]
        n_outcross = n_seeds - n_self[i]
        male_part = np.concatenate([
            np.repeat(f, n_self[i]),
            np.random.choice(males, n_outcross) if males and n_outcross > 0 else np.array([], dtype=int)
        ])
        if len(male_part) < n_seeds:
            male_part = np.concatenate([male_part, np.random.choice(males or [f], n_seeds - len(male_part))])
        female_part = np.repeat(f + 1, n_seeds)
        male_part = male_part + 1
        cross_plan.append(np.column_stack([female_part, male_part]))
    cross_plan = np.vstack(cross_plan)
    return make_cross(pop, cross_plan, 1, sim_param)


def hybrid_cross(females: Pop, males: Pop, cross_plan: str = "testcross",
                 return_hybrid_pop: bool = False,
                 sim_param: Optional[SimParam] = None) -> Pop:
    """Hybrid crossing - all combinations. From AlphaSimR hybrids.R hybridCross."""
    if sim_param is None:
        raise ValueError("simParam must be provided")
    if cross_plan == "testcross":
        cross_plan = np.column_stack([
            np.repeat(np.arange(females.n_ind), males.n_ind),
            np.tile(np.arange(males.n_ind), females.n_ind)
        ])
        if np.min(cross_plan) >= 1:
            cross_plan = cross_plan + 1
    return make_cross2(females, males, cross_plan, 1, sim_param)


def calc_gca(pop: Pop, use: str = "pheno") -> dict:
    """Calculate GCA from hybrid pop. From AlphaSimR hybrids.R calcGCA."""
    if use == "pheno":
        y = pop.pheno
    elif use == "gv":
        y = pop.gv
    elif use == "ebv":
        y = pop.ebv
    else:
        raise ValueError(f"use={use} is not a valid option")
    if y.shape[1] == 0:
        raise ValueError(f"No values for {use}")
    female = pop.mother
    male = pop.father
    female_ids = list(dict.fromkeys(female))
    male_ids = list(dict.fromkeys(male))
    n_f = len(female_ids)
    n_m = len(male_ids)
    f_idx = np.array([female_ids.index(f) for f in female])
    m_idx = np.array([male_ids.index(m) for m in male])
    y_mean = np.mean(y, axis=1)
    if y_mean.ndim == 1:
        y_mean = y_mean.reshape(-1, 1)
    n_traits = y_mean.shape[1]
    n_hyb = len(female)
    GCAf = np.zeros((n_f, n_traits))
    GCAm = np.zeros((n_m, n_traits))
    SCA = np.zeros((n_hyb, n_traits))
    for t in range(n_traits):
        grand_mean = np.mean(y_mean[:, t])
        for i in range(n_f):
            mask = f_idx == i
            GCAf[i, t] = np.mean(y_mean[mask, t]) - grand_mean
        for j in range(n_m):
            mask = m_idx == j
            GCAm[j, t] = np.mean(y_mean[mask, t]) - grand_mean
        for k in range(n_hyb):
            SCA[k, t] = y_mean[k, t] - grand_mean - GCAf[f_idx[k], t] - GCAm[m_idx[k], t]
    return {'GCAf': GCAf, 'GCAm': GCAm, 'SCA': SCA}


# Public API: snake_case module functions; SimParam methods keep AlphaSimR-style names.
__all__ = [
    "MapPop", "TraitA", "TraitAD", "TraitAE", "TraitAG", "SimParam", "Pop", "MultiPop",
    "run_macs", "quick_haplo", "new_map_pop", "import_haplo", "new_pop", "new_empty_pop", "new_multi_pop",
    "pull_seg_site_geno", "add_error", "set_pheno", "mean_g", "var_g", "mean_p",
    "select_ind", "rand_cross", "make_cross", "make_cross2", "rand_cross2", "self_", "make_dh",
    "select_cross", "select_op", "hybrid_cross", "calc_gca", "as_categorical",
    "merge_pops", "edit_genome", "import_gen_map", "get_num_threads", "sample_int",
]