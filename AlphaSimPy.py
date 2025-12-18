import numpy as np
import random
from typing import Optional, List, Union, Dict, Any
from dataclasses import dataclass
import multiprocessing as mp
from concurrent.futures import ThreadPoolExecutor
import threading

# Try to import C++ bindings, fall back to Python implementation if not available
try:
    import _alphasimpy_cpp
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

        print(command)
        
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


# Convenience function with the same name as AlphaSimR
def runMacs(nInd, nChr=1, segSites=None, inbred=False, species="GENERIC",
            split=None, ploidy=2, manualCommand=None, manualGenLen=None, nThreads=None):
    """
    Convenience function with AlphaSimR-style parameter names
    """
    return run_macs(
        n_ind=nInd,
        n_chr=nChr,
        seg_sites=segSites,
        inbred=inbred,
        species=species,
        split=split,
        ploidy=ploidy,
        manual_command=manualCommand,
        manual_gen_len=manualGenLen,
        n_threads=nThreads
    )


def run_macs2(n_ind: int,
              n_chr: int = 1,
              seg_sites: Optional[Union[int, List[int]]] = None,
              ne: float = 100,
              bp: float = 1e8,
              gen_len: float = 1,
              mut_rate: float = 2.5e-8,
              hist_ne: Optional[List[float]] = None,
              hist_gen: Optional[List[float]] = None,
              inbred: bool = False,
              split: Optional[int] = None,
              ploidy: int = 2,
              return_command: bool = False,
              n_threads: Optional[int] = None) -> Union[MapPop, str]:
    """
    Alternative wrapper for MaCS simulation
    
    A wrapper function for runMacs that provides a more intuitive interface
    for writing custom commands in MaCS. It effectively automates the creation
    of an appropriate line for the manualCommand argument in runMacs using
    user supplied variables.
    
    Parameters:
    -----------
    n_ind : int
        Number of individuals to simulate
    n_chr : int, default=1
        Number of chromosomes to simulate
    seg_sites : int or list of int, optional
        Number of segregating sites to keep per chromosome
    ne : float, default=100
        Effective population size
    bp : float, default=1e8
        Base pair length of chromosome
    gen_len : float, default=1
        Genetic length of chromosome in Morgans
    mut_rate : float, default=2.5e-8
        Per base pair mutation rate
    hist_ne : list of float, optional
        Effective population size in previous generations
    hist_gen : list of float, optional
        Number of generations ago for effective population sizes given in hist_ne
    inbred : bool, default=False
        Should founder individuals be inbred
    split : int, optional
        Optional historic population split in terms of generations ago
    ploidy : int, default=2
        Ploidy level of organism
    return_command : bool, default=False
        Should the command passed to manualCommand in runMacs be returned.
        If True, MaCS will not be called and the command is returned instead.
    n_threads : int, optional
        Number of threads for parallel simulation. If None, automatically detected.
    
    Returns:
    --------
    MapPop or str
        A MapPop object containing the simulated population, or if return_command
        is True, a string giving the MaCS command passed to the manualCommand
        argument of runMacs.
    """
    if hist_ne is None:
        hist_ne = [500, 1500, 6000, 12000, 100000]
    if hist_gen is None:
        hist_gen = [100, 1000, 10000, 100000, 1000000]
    
    if len(hist_ne) != len(hist_gen):
        raise ValueError("hist_ne and hist_gen must have the same length")
    
    # Adjust Ne according to ploidy level
    ne = ne * (ploidy / 2.0)
    if hist_ne is not None:
        hist_ne = [h * (ploidy / 2.0) for h in hist_ne]
    
    # Build species parameters
    species_params = f"{int(bp)} -t {4 * ne * mut_rate} -r {4 * ne * gen_len / bp}"
    
    # Build species history
    species_hist = ""
    if len(hist_ne) > 0:
        hist_ne_scaled = [h / ne for h in hist_ne]
        hist_gen_scaled = [g / (4 * ne) for g in hist_gen]
        for i in range(len(hist_ne_scaled)):
            species_hist += f" -eN {hist_gen_scaled[i]} {hist_ne_scaled[i]}"
    
    # Build command
    if split is None:
        command = f"{species_params}{species_hist}"
    else:
        pop_size = n_ind if inbred else ploidy * n_ind
        command = f"{species_params} -I 2 {pop_size // 2} {pop_size // 2}{species_hist} -ej {split / (4 * ne) + 0.000001} 2 1"
    
    if return_command:
        return command
    
    return run_macs(
        n_ind=n_ind,
        n_chr=n_chr,
        seg_sites=seg_sites,
        inbred=inbred,
        species="GENERIC",  # Species not used when manual_command is provided
        split=None,  # Split is handled in command
        ploidy=ploidy,
        manual_command=command,
        manual_gen_len=gen_len,
        n_threads=n_threads
    )


# Convenience function with AlphaSimR-style parameter names
def runMacs2(nInd, nChr=1, segSites=None, Ne=100, bp=1e8, genLen=1,
             mutRate=2.5e-8, histNe=None, histGen=None, inbred=False,
             split=None, ploidy=2, returnCommand=False, nThreads=None):
    """
    Convenience function with AlphaSimR-style parameter names
    """
    return run_macs2(
        n_ind=nInd,
        n_chr=nChr,
        seg_sites=segSites,
        ne=Ne,
        bp=bp,
        gen_len=genLen,
        mut_rate=mutRate,
        hist_ne=histNe,
        hist_gen=histGen,
        inbred=inbred,
        split=split,
        ploidy=ploidy,
        return_command=returnCommand,
        n_threads=nThreads
    )


# Trait classes
@dataclass
class TraitA:
    """Additive trait class"""
    qtl_loci: List[int]
    add_eff: np.ndarray
    intercept: float
    name: str


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
class TraitAD:
    """Additive + Dominance trait class"""
    qtl_loci: List[int]
    add_eff: np.ndarray
    dom_eff: np.ndarray
    intercept: float
    name: str


@dataclass
class TraitADG:
    """Additive + Dominance + GxE trait class"""
    qtl_loci: List[int]
    add_eff: np.ndarray
    dom_eff: np.ndarray
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
    
    def __init__(self, founder_pop: MapPop):
        """
        Initialize SimParam with a founder population
        
        Parameters:
        -----------
        founder_pop : MapPop
            The founder population used for variance scaling
        """
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
    
    def add_trait_adg(self, n_qtl_per_chr: Union[int, List[int]], 
                     mean: Union[float, List[float]] = 0,
                     var: Union[float, List[float]] = 1,
                     var_env: Union[float, List[float]] = 0,
                     var_gxe: Union[float, List[float]] = 1e-6,
                     mean_dd: Union[float, List[float]] = 0,
                     var_dd: Union[float, List[float]] = 0,
                     cor_a: Optional[np.ndarray] = None,
                     cor_dd: Optional[np.ndarray] = None,
                     cor_gxe: Optional[np.ndarray] = None,
                     use_var_a: bool = True,
                     gamma: Union[bool, List[bool]] = False,
                     shape: Union[float, List[float]] = 1,
                     force: bool = False,
                     name: Optional[Union[str, List[str]]] = None):
        """
        Add additive + dominance + GxE trait(s) to the simulation
        
        Parameters:
        -----------
        n_qtl_per_chr : int or list of int
            Number of QTLs per chromosome
        mean : float or list of float, default=0
            Desired mean genetic values
        var : float or list of float, default=1
            Desired genetic variances
        var_env : float or list of float, default=0
            Environmental variances
        var_gxe : float or list of float, default=1e-6
            Total genotype-by-environment variances
        mean_dd : float or list of float, default=0
            Mean dominance degree
        var_dd : float or list of float, default=0
            Variance of dominance degree
        cor_a : numpy.ndarray, optional
            Correlation matrix between additive effects
        cor_dd : numpy.ndarray, optional
            Correlation matrix between dominance degrees
        cor_gxe : numpy.ndarray, optional
            Correlation matrix between GxE effects
        use_var_a : bool, default=True
            Tune according to additive genetic variance if True
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
        if isinstance(var_env, (int, float)):
            var_env = [var_env]
        if isinstance(var_gxe, (int, float)):
            var_gxe = [var_gxe]
        if isinstance(mean_dd, (int, float)):
            mean_dd = [mean_dd]
        if isinstance(var_dd, (int, float)):
            var_dd = [var_dd]
        if isinstance(gamma, bool):
            gamma = [gamma] * len(mean)
        if isinstance(shape, (int, float)):
            shape = [shape] * len(mean)
        
        n_traits = len(mean)
        
        # Set up correlation matrices
        if cor_a is None:
            cor_a = np.eye(n_traits)
        if cor_dd is None:
            cor_dd = np.eye(n_traits)
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
        if len(mean_dd) != n_traits:
            raise ValueError("mean_dd must have same length as mean")
        if len(var_dd) != n_traits:
            raise ValueError("var_dd must have same length as mean")
        if not np.allclose(cor_a, cor_a.T):
            raise ValueError("cor_a must be symmetric")
        if not np.allclose(cor_dd, cor_dd.T):
            raise ValueError("cor_dd must be symmetric")
        if not np.allclose(cor_gxe, cor_gxe.T):
            raise ValueError("cor_gxe must be symmetric")
        if cor_a.shape[0] != n_traits:
            raise ValueError("cor_a must be square with dimension equal to number of traits")
        if cor_dd.shape[0] != n_traits:
            raise ValueError("cor_dd must be square with dimension equal to number of traits")
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
        
        # Sample dominance effects
        dom_eff = self._samp_dom_eff(qtl_loci, n_traits, add_eff, cor_dd, mean_dd, var_dd)
        
        # Sample GxE effects
        gxe_eff = self._samp_add_eff(qtl_loci, n_traits, cor_gxe, gamma=[False] * n_traits, shape=shape)
        
        # Create traits
        for i in range(n_traits):
            # Create base TraitAD
            trait = TraitAD(
                qtl_loci=flat_qtl_loci,
                add_eff=add_eff[:, i],
                dom_eff=dom_eff[:, i],
                intercept=0,
                name=name[i]
            )
            
            # Calculate genetic parameters
            tmp = self._calc_gen_param_ad(trait, self.founder_pop)
            
            if use_var_a:
                bv_var = self._pop_var(tmp['bv'])
                if np.isscalar(bv_var):
                    scale = np.sqrt(var[i]) / np.sqrt(bv_var)
                else:
                    scale = np.sqrt(var[i]) / np.sqrt(bv_var[0])
            else:
                gv_var = self._pop_var(tmp['gv'])
                if np.isscalar(gv_var):
                    scale = np.sqrt(var[i]) / np.sqrt(gv_var)
                else:
                    scale = np.sqrt(var[i]) / np.sqrt(gv_var[0])
            
            trait.add_eff = trait.add_eff * scale
            trait.dom_eff = trait.dom_eff * scale
            trait.intercept = mean[i] - np.mean(tmp['gv'] * scale)
            
            # GxE component
            trait_g = TraitA(
                qtl_loci=flat_qtl_loci,
                add_eff=gxe_eff[:, i],
                intercept=0,
                name=name[i]
            )
            tmp_g = self._calc_gen_param(trait_g, self.founder_pop)
            
            gv_var_g = self._pop_var(tmp_g['gv'])
            if np.isscalar(gv_var_g):
                gv_var_g_val = gv_var_g
            else:
                gv_var_g_val = gv_var_g[0]
            
            if var_env[i] == 0:
                scale_g = np.sqrt(var_gxe[i]) / np.sqrt(gv_var_g_val)
                trait_adg = TraitADG(
                    qtl_loci=flat_qtl_loci,
                    add_eff=trait.add_eff,
                    dom_eff=trait.dom_eff,
                    intercept=trait.intercept,
                    name=name[i],
                    gxe_eff=gxe_eff[:, i] * scale_g,
                    gxe_int=0 - np.mean(tmp_g['gv'] * scale_g),
                    env_var=1
                )
            else:
                scale_g = np.sqrt(var_gxe[i] / var_env[i]) / np.sqrt(gv_var_g_val)
                trait_adg = TraitADG(
                    qtl_loci=flat_qtl_loci,
                    add_eff=trait.add_eff,
                    dom_eff=trait.dom_eff,
                    intercept=trait.intercept,
                    name=name[i],
                    gxe_eff=gxe_eff[:, i] * scale_g,
                    gxe_int=1 - np.mean(tmp_g['gv'] * scale_g),
                    env_var=var_env[i]
                )
            
            if use_var_a:
                gv_var_scaled = self._pop_var(tmp['gv'] * scale)
                if np.isscalar(gv_var_scaled):
                    self._add_trait(trait_adg, var[i], gv_var_scaled)
                else:
                    self._add_trait(trait_adg, var[i], gv_var_scaled[0])
            else:
                bv_var_scaled = self._pop_var(tmp['bv'] * scale)
                if np.isscalar(bv_var_scaled):
                    self._add_trait(trait_adg, bv_var_scaled, var[i])
                else:
                    self._add_trait(trait_adg, bv_var_scaled[0], var[i])
        
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
    
    def _samp_dom_eff(self, qtl_loci: List[List[int]], n_traits: int,
                     add_eff: np.ndarray, cor_dd: np.ndarray,
                     mean_dd: List[float], var_dd: List[float]) -> np.ndarray:
        """
        Sample dominance effects for QTL
        
        Parameters:
        -----------
        qtl_loci : list of list of int
            QTL loci for each chromosome
        n_traits : int
            Number of traits
        add_eff : numpy.ndarray
            Previously sampled additive effects
        cor_dd : numpy.ndarray
            Correlation matrix between dominance degrees
        mean_dd : list of float
            Mean value of dominance degrees
        var_dd : list of float
            Variance of dominance degrees
        
        Returns:
        --------
        numpy.ndarray
            Dominance effects matrix
        """
        n_qtl = sum(len(chr_loci) for chr_loci in qtl_loci)
        
        # Sample dominance degrees from normal distribution with correlation
        chol = np.linalg.cholesky(cor_dd)
        dom_eff = np.random.normal(0, 1, (n_qtl, n_traits))
        dom_eff = dom_eff @ chol.T
        
        # Scale by variance and add mean
        for i in range(n_traits):
            dom_eff[:, i] = dom_eff[:, i] * np.sqrt(var_dd[i]) + mean_dd[i]
        
        # Multiply by absolute additive effects
        dom_eff = np.abs(add_eff) * dom_eff
        
        return dom_eff
    
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
    
    def _calc_gen_param_ad(self, trait: TraitAD, pop: MapPop) -> Dict[str, np.ndarray]:
        """
        Calculate genetic parameters for a trait with additive and dominance effects
        
        Parameters:
        -----------
        trait : TraitAD
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
        bv = np.zeros(n_ind)  # Breeding values (additive only)
        
        # Calculate genetic values for each individual
        for ind_idx in range(n_ind):
            total_gv = trait.intercept
            total_bv = 0  # Breeding value (additive only)
            
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
                    
                    # Add QTL additive effect
                    if qtl_idx < len(trait.add_eff):
                        add_contribution = trait.add_eff[qtl_idx] * scaled_dosage
                        total_gv += add_contribution
                        total_bv += add_contribution
                    
                    # Add dominance effect (only for heterozygous genotypes)
                    if qtl_idx < len(trait.dom_eff):
                        # For diploid: dominance only when dosage = 1 (heterozygous)
                        # For polyploid: dominance when not fully homozygous
                        if pop.ploidy == 2:
                            if dosage == 1:  # Heterozygous
                                total_gv += trait.dom_eff[qtl_idx]
                        else:
                            # For polyploids, use a simplified dominance model
                            # Dominance effect is proportional to heterozygosity
                            heterozygosity = 1.0 - abs(scaled_dosage)  # 1 when heterozygous, 0 when homozygous
                            total_gv += trait.dom_eff[qtl_idx] * heterozygosity
            
            gv[ind_idx] = total_gv
            bv[ind_idx] = total_bv
        
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
    
    def set_var_e(self, h2=None, H2=None, var_e=None, cor_e=None):
        """
        Set error variance for traits
        
        Parameters:
        -----------
        h2 : float or list of float, optional
            Narrow-sense heritabilities for each trait
        H2 : float or list of float, optional
            Broad-sense heritabilities for each trait
        var_e : float, list of float, or numpy.ndarray, optional
            Error variance(s) for each trait. Can be a vector or covariance matrix.
        cor_e : numpy.ndarray, optional
            Correlation matrix for error variances
        
        Returns:
        --------
        SimParam
            Returns self for method chaining
        """
        n_traits = self.n_traits
        
        # Check validity of corE, if supplied
        if cor_e is not None:
            if not np.allclose(cor_e, cor_e.T):
                raise ValueError("corE must be symmetric")
            if cor_e.shape[0] != n_traits:
                raise ValueError(f"corE must be square with dimension equal to number of traits ({n_traits})")
        
        # Set error variances
        if h2 is not None:
            if isinstance(h2, (int, float)):
                h2 = [h2]
            if len(h2) != n_traits:
                raise ValueError(f"h2 must have length equal to number of traits ({n_traits})")
            if not all(self._var_g[i] > 0 for i in range(n_traits)):
                raise ValueError("All varG must be > 0 when using h2")
            if not all(self._var_a[i] > 0 for i in range(n_traits)):
                raise ValueError("All varA must be > 0 when using h2")
            
            var_e = []
            for i in range(n_traits):
                tmp = self._var_a[i] / h2[i] - self._var_g[i]
                if tmp < 0:
                    raise ValueError(f"h2={h2[i]} is not possible for trait {i+1}")
                var_e.append(tmp)
            self._var_e = var_e
            
        elif H2 is not None:
            if isinstance(H2, (int, float)):
                H2 = [H2]
            if len(H2) != n_traits:
                raise ValueError(f"H2 must have length equal to number of traits ({n_traits})")
            
            var_e = []
            for i in range(n_traits):
                tmp = self._var_g[i] / H2[i] - self._var_g[i]
                var_e.append(tmp)
            self._var_e = var_e
            
        elif var_e is not None:
            if isinstance(var_e, (int, float)):
                var_e = [var_e] * n_traits
            elif isinstance(var_e, np.ndarray):
                if var_e.ndim == 2:
                    # Matrix - check dimensions
                    if var_e.shape[0] != n_traits or var_e.shape[1] != n_traits:
                        raise ValueError(f"varE matrix must be {n_traits}x{n_traits}")
                    self._var_e = var_e
                else:
                    # Vector
                    if len(var_e) != n_traits:
                        raise ValueError(f"varE must have length equal to number of traits ({n_traits})")
                    self._var_e = var_e.tolist()
            else:
                # List
                if len(var_e) != n_traits:
                    raise ValueError(f"varE must have length equal to number of traits ({n_traits})")
                self._var_e = var_e
        else:
            self._var_e = [np.nan] * n_traits
        
        # Set error correlations
        if cor_e is not None:
            if isinstance(self._var_e, np.ndarray) and self._var_e.ndim == 2:
                var_e_diag = np.diag(self._var_e)
            else:
                var_e_diag = np.array(self._var_e)
            
            var_e_sqrt = np.diag(np.sqrt(var_e_diag))
            var_e_matrix = var_e_sqrt @ cor_e @ var_e_sqrt
            self._var_e = var_e_matrix
        
        return self


# Add convenience methods to SimParam class
def addTraitA(self, nQtlPerChr, mean=0, var=1, corA=None, gamma=False, shape=1, force=False, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_trait_a(nQtlPerChr, mean, var, corA, gamma, shape, force, name)

def addTraitAG(self, nQtlPerChr, mean=0, var=1, varGxE=1e-6, varEnv=0, corA=None, corGxE=None, gamma=False, shape=1, force=False, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_trait_ag(nQtlPerChr, mean, var, varGxE, varEnv, corA, corGxE, gamma, shape, force, name)

def addTraitADG(self, nQtlPerChr, mean=0, var=1, varEnv=0, varGxE=1e-6, meanDD=0, varDD=0, corA=None, corDD=None, corGxE=None, useVarA=True, gamma=False, shape=1, force=False, name=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.add_trait_adg(nQtlPerChr, mean, var, varEnv, varGxE, meanDD, varDD, corA, corDD, corGxE, useVarA, gamma, shape, force, name)

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

def setVarE(self, h2=None, H2=None, varE=None, corE=None):
    """Convenience method with AlphaSimR-style parameter names"""
    return self.set_var_e(h2, H2, varE, corE)

# Add methods to SimParam class
SimParam.addTraitA = addTraitA
SimParam.addTraitAG = addTraitAG
SimParam.addTraitADG = addTraitADG
SimParam.restrSegSites = restrSegSites
SimParam.addSnpChip = addSnpChip
SimParam.setTrackPed = setTrackPed
SimParam.setTrackRec = setTrackRec
SimParam.resetPed = resetPed
SimParam.setSexes = setSexes
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


@dataclass
class HybridPop:
    """Hybrid population class - lightweight version of Pop without genotypic data"""
    n_ind: int
    id: List[str]
    mother: List[str]
    father: List[str]
    n_traits: int
    gv: np.ndarray
    pheno: np.ndarray
    gxe: List[Any]
    
    def __post_init__(self):
        """Validate the HybridPop object"""
        if any(' ' in individual_id for individual_id in self.id):
            raise ValueError("id cannot contain spaces")
        if any(' ' in mother_id for mother_id in self.mother):
            raise ValueError("mother cannot contain spaces")
        if any(' ' in father_id for father_id in self.father):
            raise ValueError("father cannot contain spaces")
        if self.n_ind != len(self.id):
            raise ValueError("nInd != length(id)")
        if self.n_ind != len(self.mother):
            raise ValueError("nInd != length(mother)")
        if self.n_ind != len(self.father):
            raise ValueError("nInd != length(father)")
        if self.n_ind != self.gv.shape[0]:
            raise ValueError("nInd != nrow(gv)")
        if self.n_ind != self.pheno.shape[0]:
            raise ValueError("nInd != nrow(pheno)")
        if self.n_traits != self.gv.shape[1]:
            raise ValueError("nTraits != ncol(gv)")
        if self.n_traits != self.pheno.shape[1]:
            raise ValueError("nTraits != ncol(pheno)")
        if self.n_traits != len(self.gxe):
            raise ValueError("nTraits != length(gxe)")


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


# Convenience function with AlphaSimR-style name
def newPop(raw_pop: MapPop, sim_param: Optional[SimParam] = None, **kwargs) -> Pop:
    """
    Convenience function with AlphaSimR-style parameter names
    """
    return new_pop(raw_pop, sim_param, **kwargs)


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


# Convenience functions with AlphaSimR-style names
def meanG(pop: Pop) -> np.ndarray:
    """Convenience function with AlphaSimR-style name"""
    return mean_g(pop)


def varG(pop: Pop) -> np.ndarray:
    """Convenience function with AlphaSimR-style name"""
    return var_g(pop)


def meanP(pop: Pop) -> np.ndarray:
    """Convenience function with AlphaSimR-style name"""
    return mean_p(pop)


def nInd(pop: Pop) -> int:
    """
    Return the number of individuals in a population.

    Mirrors AlphaSimR::nInd by exposing the n_ind field.
    """
    return int(pop.n_ind)


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


def setPheno(pop: Pop, varE: Optional[Union[float, List[float], np.ndarray]] = None,
             reps: int = 1, simParam: Optional[SimParam] = None) -> Pop:
    """Convenience function with AlphaSimR-style parameter names"""
    return set_pheno(pop, varE, reps, simParam)


def merge_pops(pop_list: List[Pop]) -> Pop:
    """
    Merge a list of Pop objects into a single Pop.

    This function closely mirrors AlphaSimR::mergePops for single Pop
    objects and assumes all populations share the same genetic map.
    """
    # Drop None entries to mirror R's behavior of skipping NULL
    pop_list = [p for p in pop_list if p is not None]
    if not pop_list:
        raise ValueError("popList must contain at least one Pop")

    # All elements must be Pop instances
    for p in pop_list:
        if not isinstance(p, Pop):
            raise TypeError("All elements of popList must be Pop instances")

    first = pop_list[0]

    # Check core structural compatibility: n_chr, ploidy, n_loci
    for p in pop_list[1:]:
        if p.n_chr != first.n_chr:
            raise ValueError("All populations must have the same n_chr")
        if p.ploidy != first.ploidy:
            raise ValueError("All populations must have the same ploidy")
        if p.n_loci != first.n_loci:
            raise ValueError("All populations must have identical n_loci")

    # Check genetic map compatibility (not present in AlphaSimR::mergePops but
    # required for the Python Pop definition)
    for p in pop_list[1:]:
        for chr_idx in range(first.n_chr):
            if not np.array_equal(first.gen_map[chr_idx], p.gen_map[chr_idx]):
                raise ValueError("All populations must have identical gen_map")
            if first.centromere[chr_idx] != p.centromere[chr_idx]:
                raise ValueError("All populations must have identical centromere positions")

    # n_traits
    n_traits_list = [p.n_traits for p in pop_list]
    if any(n != n_traits_list[0] for n in n_traits_list):
        raise ValueError("All populations must have the same n_traits")
    n_traits = n_traits_list[0]

    # n_ind per population and total
    n_ind_list = [p.n_ind for p in pop_list]
    total_n_ind = int(sum(n_ind_list))

    # id, iid, parents, sex, fixed effects
    id_list: List[str] = []
    iid_list: List[int] = []
    mother_list: List[str] = []
    father_list: List[str] = []
    sex_list: List[str] = []
    fix_eff_list: List[int] = []

    for p in pop_list:
        id_list.extend(p.id)
        iid_list.extend(p.iid)
        mother_list.extend(p.mother)
        father_list.extend(p.father)
        sex_list.extend(p.sex)
        fix_eff_list.extend(p.fix_eff)

    # misc: only merge when all populations have the same non-empty key set
    misc_merged: Dict[str, Any] = {}
    if all(isinstance(p.misc, dict) for p in pop_list):
        key_sets = [set(p.misc.keys()) for p in pop_list]
        if key_sets and all(ks == key_sets[0] for ks in key_sets) and len(key_sets[0]) > 0:
            for key in key_sets[0]:
                vals = [p.misc[key] for p in pop_list]
                first_val = vals[0]
                # If matrix-like (2D ndarray), stack row-wise as in rbind
                if isinstance(first_val, np.ndarray) and first_val.ndim == 2:
                    misc_merged[key] = np.vstack(vals)
                else:
                    # Otherwise, concatenate sensibly
                    try:
                        misc_merged[key] = np.concatenate(vals)
                    except Exception:
                        combined: List[Any] = []
                        for v in vals:
                            if isinstance(v, (list, tuple)):
                                combined.extend(v)
                            else:
                                combined.append(v)
                        misc_merged[key] = combined

    # n_traits, gv, pheno
    if n_traits > 0 and total_n_ind > 0:
        gv = np.vstack([p.gv for p in pop_list])
        pheno = np.vstack([p.pheno for p in pop_list])
    else:
        gv = np.empty((total_n_ind, n_traits))
        pheno = np.empty((total_n_ind, n_traits))

    # ebv: only merge if all have same number of columns, otherwise 0-col matrix
    ebv_cols = [p.ebv.shape[1] for p in pop_list]
    if all(c == ebv_cols[0] for c in ebv_cols):
        if ebv_cols[0] > 0 and total_n_ind > 0:
            ebv = np.vstack([p.ebv for p in pop_list])
        else:
            ebv = np.empty((total_n_ind, 0))
    else:
        ebv = np.empty((total_n_ind, 0))

    # gxe: concatenate per trait when present
    if n_traits >= 1:
        gxe_merged: List[Any] = [None] * n_traits
        for trait_idx in range(n_traits):
            first_val = pop_list[0].gxe[trait_idx] if pop_list[0].gxe else None
            if first_val is None:
                continue
            vals = [p.gxe[trait_idx] for p in pop_list]
            if isinstance(first_val, np.ndarray):
                gxe_merged[trait_idx] = np.concatenate(vals, axis=0)
            else:
                combined: List[Any] = []
                for v in vals:
                    if isinstance(v, (list, tuple)):
                        combined.extend(v)
                    else:
                        combined.append(v)
                gxe_merged[trait_idx] = combined
    else:
        gxe_merged = []

    # geno: concatenate along individual dimension for each chromosome
    merged_geno: List[np.ndarray] = []
    for chr_idx in range(first.n_chr):
        geno_arrays = [p.geno[chr_idx] for p in pop_list]
        base_shape = geno_arrays[0].shape[:2]
        for arr in geno_arrays:
            if arr.shape[:2] != base_shape:
                raise ValueError("Incompatible geno shapes across populations")
        if total_n_ind > 0:
            merged_chr = np.concatenate(geno_arrays, axis=2)
        else:
            merged_chr = geno_arrays[0][:, :, :0]
        merged_geno.append(merged_chr)

    # inbred flag: AND across populations (as in AlphaSimR for MapPop)
    inbred_flag = all(p.inbred for p in pop_list)

    # Construct merged population
    merged_pop = Pop(
        n_ind=total_n_ind,
        n_chr=first.n_chr,
        ploidy=first.ploidy,
        n_loci=first.n_loci,
        geno=merged_geno,
        gen_map=first.gen_map,
        centromere=first.centromere,
        inbred=inbred_flag,
        id=id_list,
        iid=iid_list,
        mother=mother_list,
        father=father_list,
        sex=sex_list,
        n_traits=n_traits,
        gv=gv,
        pheno=pheno,
        ebv=ebv,
        gxe=gxe_merged,
        fix_eff=fix_eff_list,
        misc=misc_merged,
        misc_pop={}
    )

    return merged_pop


def mergePops(popList: List[Pop]) -> Pop:
    """Convenience wrapper with AlphaSimR-style name."""
    return merge_pops(popList)


# Selection functions
def select_ind(pop: Pop, n_ind: int, trait: Union[int, callable] = 1,
               use: Union[str, callable] = "pheno", sex: str = "B",
               select_top: bool = True, return_pop: bool = True,
               candidates: Optional[List[int]] = None,
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
    
    # Select individuals
    eligible_values = values[eligible]
    sorted_indices = np.argsort(eligible_values)
    if select_top:
        selected_indices = sorted_indices[-n_ind:]
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


def selectInd(pop: Pop, nInd: int, trait: Union[int, callable] = 1,
              use: str = "pheno", sex: str = "B", selectTop: bool = True,
              returnPop: bool = True, candidates: Optional[List[int]] = None,
              simParam: Optional[SimParam] = None, **kwargs) -> Union[Pop, List[int]]:
    """Convenience function with AlphaSimR-style parameter names"""
    return select_ind(pop, nInd, trait, use, sex, selectTop, returnPop, candidates, simParam, **kwargs)


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


def randCross(pop: Pop, nCrosses: int, nProgeny: int = 1, balance: bool = True,
              parents: Optional[List[int]] = None, ignoreSexes: bool = False,
              simParam: Optional[SimParam] = None) -> Pop:
    """Convenience function with AlphaSimR-style parameter names"""
    return rand_cross(pop, nCrosses, nProgeny, balance, parents, ignoreSexes, simParam)


def make_cross2(females: Pop, males: Pop, cross_plan: np.ndarray,
                sim_param: Optional[SimParam] = None) -> Pop:
    """
    Make crosses between two populations using a cross plan
    
    Parameters:
    -----------
    females : Pop
        Female population
    males : Pop
        Male population
    cross_plan : numpy.ndarray
        Matrix with two columns representing female and male parent indices
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop
        New population of progeny
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    if females.ploidy % 2 != 0 or males.ploidy % 2 != 0:
        raise ValueError("You can not cross individuals with odd ploidy levels")
    
    n_crosses = cross_plan.shape[0]
    if cross_plan.shape[1] != 2:
        raise ValueError("crossPlan must have 2 columns")
    
    # Validate cross plan indices
    if np.max(cross_plan[:, 0]) >= females.n_ind or np.min(cross_plan[:, 0]) < 0:
        raise ValueError("Invalid crossPlan: female indices out of range")
    if np.max(cross_plan[:, 1]) >= males.n_ind or np.min(cross_plan[:, 1]) < 0:
        raise ValueError("Invalid crossPlan: male indices out of range")
    
    # Generate progeny genotypes
    total_progeny = n_crosses
    progeny_geno = []
    
    for chr_idx in range(females.n_chr):
        n_packed = females.geno[chr_idx].shape[0]
        chr_geno = np.zeros((n_packed, females.ploidy, total_progeny), dtype=np.uint8)
        
        for cross_idx in range(n_crosses):
            female_idx = int(cross_plan[cross_idx, 0])
            male_idx = int(cross_plan[cross_idx, 1])
            
            # Simple crossing simulation - randomly inherit from each parent
            for packed_idx in range(n_packed):
                for haplo in range(females.ploidy):
                    if np.random.random() < 0.5:
                        # Inherit from female
                        chr_geno[packed_idx, haplo, cross_idx] = females.geno[chr_idx][packed_idx, haplo, female_idx]
                    else:
                        # Inherit from male
                        chr_geno[packed_idx, haplo, cross_idx] = males.geno[chr_idx][packed_idx, haplo, male_idx]
        
        progeny_geno.append(chr_geno)
    
    # Create progeny population
    female_parents = [females.id[int(cross_plan[i, 0])] for i in range(n_crosses)]
    male_parents = [males.id[int(cross_plan[i, 1])] for i in range(n_crosses)]
    
    progeny_pop = Pop(
        n_ind=total_progeny, n_chr=females.n_chr, ploidy=females.ploidy, n_loci=females.n_loci,
        geno=progeny_geno, gen_map=females.gen_map, centromere=females.centromere,
        inbred=False,  # Hybrids are not inbred
        id=[f"{female_parents[i]}_{male_parents[i]}" for i in range(total_progeny)],
        iid=list(range(sim_param._last_id + 1, sim_param._last_id + total_progeny + 1)),
        mother=female_parents,
        father=male_parents,
        sex=["H"] * total_progeny,
        n_traits=females.n_traits,
        gv=np.zeros((total_progeny, females.n_traits)),
        pheno=np.zeros((total_progeny, females.n_traits)),
        ebv=np.empty((total_progeny, 0)),
        gxe=[None] * females.n_traits,
        fix_eff=[1] * total_progeny,
        misc={},
        misc_pop={}
    )
    
    # Calculate genetic values for progeny
    if females.n_traits > 0:
        progeny_pop.gv = _get_gv_index(progeny_pop, sim_param)
        progeny_pop.pheno = progeny_pop.gv.copy()
    
    # Update last ID
    sim_param._last_id += total_progeny
    
    return progeny_pop


def calc_coef(X: np.ndarray, Y: np.ndarray) -> np.ndarray:
    """
    Calculate coefficients by solving X * coef = Y
    Equivalent to AlphaSimR's calcCoef function
    
    Parameters:
    -----------
    X : numpy.ndarray
        Design matrix
    Y : numpy.ndarray
        Response matrix
    
    Returns:
    --------
    numpy.ndarray
        Coefficient matrix
    """
    # Use least squares if system is overdetermined or singular
    try:
        return np.linalg.solve(X, Y)
    except np.linalg.LinAlgError:
        # Use least squares if solve fails
        return np.linalg.lstsq(X, Y, rcond=None)[0]


def _get_hybrid_gv(females: Pop, female_parents: np.ndarray,
                   males: Pop, male_parents: np.ndarray,
                   sim_param: SimParam) -> tuple:
    """
    Calculate hybrid genetic values
    
    Parameters:
    -----------
    females : Pop
        Female population
    female_parents : numpy.ndarray
        Indices of female parents
    males : Pop
        Male population
    male_parents : numpy.ndarray
        Indices of male parents
    sim_param : SimParam
        Simulation parameters
    
    Returns:
    --------
    tuple
        (gv, gxe) where gv is genetic values and gxe is GxE slopes (if applicable)
    """
    n_hybrids = len(female_parents)
    n_traits = sim_param.n_traits
    
    if n_traits == 0:
        return (np.array([]), [])
    
    # Create temporary hybrid population once for all traits
    # In a full implementation, this would use the C++ getHybridGv function
    # For now, we'll approximate by creating actual crosses
    temp_cross_plan = np.column_stack([female_parents, male_parents])
    temp_hybrids = make_cross2(females, males, temp_cross_plan, sim_param)
    
    gv = temp_hybrids.gv
    gxe = [None] * n_traits
    
    return (gv, gxe)


def hybrid_cross(females: Pop, males: Pop,
                cross_plan: Union[str, np.ndarray] = "testcross",
                return_hybrid_pop: bool = False,
                sim_param: Optional[SimParam] = None) -> Union[Pop, HybridPop]:
    """
    Hybrid crossing function for plant breeding simulations
    
    Parameters:
    -----------
    females : Pop
        Female population
    males : Pop
        Male population
    cross_plan : str or numpy.ndarray, default="testcross"
        Either "testcross" for all possible combinations or a matrix with two columns
    return_hybrid_pop : bool, default=False
        Should results be returned as HybridPop. If False returns as Pop.
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop or HybridPop
        New population of hybrids
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    if females.ploidy % 2 != 0 or males.ploidy % 2 != 0:
        raise ValueError("You can not cross individuals with odd ploidy levels")
    
    # Handle cross plan
    if isinstance(cross_plan, str):
        if cross_plan == "testcross":
            # Create all possible combinations
            cross_plan = np.array([[i, j] for i in range(females.n_ind) 
                                   for j in range(males.n_ind)])
        else:
            raise ValueError(f"crossPlan={cross_plan} is not a valid option")
    elif isinstance(cross_plan, np.ndarray):
        # Use provided cross plan
        pass
    else:
        raise ValueError("crossPlan must be 'testcross' or a numpy array")
    
    # Set IDs
    female_parents = [females.id[int(cross_plan[i, 0])] for i in range(cross_plan.shape[0])]
    male_parents = [males.id[int(cross_plan[i, 1])] for i in range(cross_plan.shape[0])]
    id_list = [f"{female_parents[i]}_{male_parents[i]}" for i in range(len(female_parents))]
    
    # Return Pop-class
    if not return_hybrid_pop:
        return make_cross2(females, males, cross_plan, sim_param)
    
    # Return HybridPop-class
    n_hybrids = len(id_list)
    n_traits = sim_param.n_traits
    
    gv = np.zeros((n_hybrids, n_traits))
    gxe = [None] * n_traits
    
    # Calculate hybrid genetic values
    if n_traits > 0:
        female_parents_idx = cross_plan[:, 0].astype(int)
        male_parents_idx = cross_plan[:, 1].astype(int)
        gv_result, gxe_result = _get_hybrid_gv(females, female_parents_idx,
                                               males, male_parents_idx,
                                               sim_param)
        gv = gv_result
        gxe = gxe_result
    
    # Add error to get phenotypes
    if n_traits > 0:
        var_e = sim_param._var_e
        if isinstance(var_e, np.ndarray) and var_e.ndim == 2:
            # Covariance matrix
            error = np.random.multivariate_normal(np.zeros(n_traits), var_e, n_hybrids)
        else:
            # Vector of variances
            error = np.random.normal(0, np.sqrt(var_e), (n_hybrids, n_traits))
        pheno = gv + error
    else:
        pheno = gv
    
    output = HybridPop(
        n_ind=n_hybrids,
        id=id_list,
        mother=female_parents,
        father=male_parents,
        n_traits=n_traits,
        gv=gv,
        pheno=pheno,
        gxe=gxe
    )
    return output


def calc_gca(pop: Union[Pop, HybridPop], use: str = "pheno") -> Dict[str, np.ndarray]:
    """
    Calculate general combining ability of test crosses
    
    Parameters:
    -----------
    pop : Pop or HybridPop
        Population with hybrid crosses
    use : str, default="pheno"
        Tabulate either genetic values "gv", estimated breeding values "ebv", or phenotypes "pheno"
    
    Returns:
    --------
    dict
        Dictionary with keys 'GCAf', 'GCAm', 'SCA' containing GCA for females, males, and SCA
    """
    use = use.lower()
    
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
    
    # Convert to factors
    female = pop.mother
    male = pop.father
    
    unique_females = list(dict.fromkeys(female))  # Preserve order
    unique_males = list(dict.fromkeys(male))
    
    female_factor = [unique_females.index(f) for f in female]
    male_factor = [unique_males.index(m) for m in male]
    
    n_females = len(unique_females)
    n_males = len(unique_males)
    n_traits = y.shape[1]
    
    # Check for balance
    if n_females == 1 or n_males == 1:
        balanced = True
    else:
        # Check if balanced design
        counts = {}
        for f, m in zip(female_factor, male_factor):
            counts[(f, m)] = counts.get((f, m), 0) + 1
        balanced = len(set(counts.values())) == 1
    
    # Calculate SCA combinations
    sca_combos = [f"{f}_{m}" for f, m in zip(female, male)]
    unique_sca = list(dict.fromkeys(sca_combos))
    sca_factor = [unique_sca.index(s) for s in sca_combos]
    n_sca = len(unique_sca)
    
    # Female GCA
    if n_females == 1:
        GCAf = np.tile(np.mean(y, axis=0), (1, 1))
    else:
        if n_males == 1:
            GCAf = y
        else:
            if balanced:
                # Calculate simple means
                GCAf = np.zeros((n_females, n_traits))
                for i, f in enumerate(unique_females):
                    mask = np.array([j == i for j in female_factor])
                    GCAf[i, :] = np.mean(y[mask, :], axis=0)
            else:
                # Calculate population marginal means using linear model
                # Create design matrix: ~female+male-1 with sum contrasts for male
                # In sum contrasts, the last level is represented as -1 in all other columns
                X = np.zeros((len(female), n_females + n_males - 1))
                for i, (f_idx, m_idx) in enumerate(zip(female_factor, male_factor)):
                    X[i, f_idx] = 1
                    if m_idx < n_males - 1:
                        X[i, n_females + m_idx] = 1
                    else:
                        # Last male level: set all male columns to -1
                        X[i, n_females:n_females + n_males - 1] = -1
                
                coef = calc_coef(X.T @ X, X.T @ y)
                GCAf = coef[:n_females, :]
    
    GCAf_df = np.column_stack([np.array(unique_females).reshape(-1, 1), GCAf])
    
    # Male GCA
    if n_males == 1:
        GCAm = np.tile(np.mean(y, axis=0), (1, 1))
    else:
        if n_females == 1:
            GCAm = y
        else:
            if balanced:
                # Calculate simple means
                GCAm = np.zeros((n_males, n_traits))
                for i, m in enumerate(unique_males):
                    mask = np.array([j == i for j in male_factor])
                    GCAm[i, :] = np.mean(y[mask, :], axis=0)
            else:
                # Calculate population marginal means using linear model
                # Create design matrix: ~male+female-1 with sum contrasts for female
                # In sum contrasts, the last level is represented as -1 in all other columns
                X = np.zeros((len(male), n_males + n_females - 1))
                for i, (m_idx, f_idx) in enumerate(zip(male_factor, female_factor)):
                    X[i, m_idx] = 1
                    if f_idx < n_females - 1:
                        X[i, n_males + f_idx] = 1
                    else:
                        # Last female level: set all female columns to -1
                        X[i, n_males:n_males + n_females - 1] = -1
                
                coef = calc_coef(X.T @ X, X.T @ y)
                GCAm = coef[:n_males, :]
    
    GCAm_df = np.column_stack([np.array(unique_males).reshape(-1, 1), GCAm])
    
    # SCA
    if n_sca == pop.n_ind:
        SCA = y
    else:
        # Calculate simple means
        SCA = np.zeros((n_sca, n_traits))
        for i, sca_combo in enumerate(unique_sca):
            mask = np.array([j == i for j in sca_factor])
            SCA[i, :] = np.mean(y[mask, :], axis=0)
    
    SCA_df = np.column_stack([np.array(unique_sca).reshape(-1, 1), SCA])
    
    return {
        'GCAf': GCAf_df,
        'GCAm': GCAm_df,
        'SCA': SCA_df
    }


def set_pheno_gca(pop: Pop, testers: Pop, use: str = "pheno",
                  h2: Optional[Union[float, List[float]]] = None,
                  H2: Optional[Union[float, List[float]]] = None,
                  var_e: Optional[Union[float, List[float], np.ndarray]] = None,
                  cor_e: Optional[np.ndarray] = None,
                  reps: int = 1, fix_eff: int = 1, p: Optional[float] = None,
                  inbred: bool = False, only_pheno: bool = False,
                  sim_param: Optional[SimParam] = None) -> Union[Pop, np.ndarray]:
    """
    Set phenotypes using general combining ability
    
    Parameters:
    -----------
    pop : Pop
        Population to set phenotypes for
    testers : Pop
        Tester population
    use : str, default="pheno"
        True genetic value ("gv") or phenotypes ("pheno")
    h2 : float or list of float, optional
        Narrow-sense heritabilities
    H2 : float or list of float, optional
        Broad-sense heritabilities
    var_e : float, list of float, or numpy.ndarray, optional
        Error variances
    cor_e : numpy.ndarray, optional
        Error correlation matrix
    reps : int, default=1
        Number of replications
    fix_eff : int, default=1
        Fixed effect to assign
    p : float, optional
        p-value for environmental covariate (GxE traits)
    inbred : bool, default=False
        Are both pop and testers fully inbred
    only_pheno : bool, default=False
        Should only phenotype be returned
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop or numpy.ndarray
        Population with updated phenotypes or phenotype matrix if only_pheno=True
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    if any(pop.id[i] in pop.id[i+1:] for i in range(len(pop.id))):
        raise ValueError("This function does not work with duplicate IDs")
    
    use = use.lower()
    
    # Make hybrids
    tmp = hybrid_cross(females=pop, males=testers, cross_plan="testcross",
                      return_hybrid_pop=inbred, sim_param=sim_param)
    
    # Get response
    if use == "pheno":
        if isinstance(tmp, HybridPop):
            # HybridPop already has phenotypes calculated
            y = tmp.pheno
        else:
            # For Pop, need to set phenotypes
            # Note: set_pheno doesn't currently support h2/H2/corE, so we'll use var_e if provided
            if var_e is not None:
                tmp = set_pheno(tmp, var_e=var_e, reps=reps, sim_param=sim_param)
            else:
                tmp = set_pheno(tmp, reps=reps, sim_param=sim_param)
            y = tmp.pheno
    elif use == "gv":
        y = tmp.gv
    else:
        raise ValueError(f"use={use} is not a valid option")
    
    if y.shape[1] == 0:
        raise ValueError(f"No values for {use}")
    
    # Calculate GCA for females
    female = tmp.mother
    unique_females = list(dict.fromkeys(female))
    
    if len(unique_females) == 1:
        GCAf = np.tile(np.mean(y, axis=0), (1, 1))
    else:
        if testers.n_ind == 1:
            GCAf = y
        else:
            # Calculate simple means
            female_factor = [unique_females.index(f) for f in female]
            n_females = len(unique_females)
            GCAf = np.zeros((n_females, y.shape[1]))
            for i, f in enumerate(unique_females):
                mask = np.array([j == i for j in female_factor])
                GCAf[i, :] = np.mean(y[mask, :], axis=0)
    
    if only_pheno:
        return GCAf
    
    pop.pheno = GCAf
    pop.fix_eff = [fix_eff] * pop.n_ind
    return pop


# Convenience functions with AlphaSimR-style names
def hybridCross(females: Pop, males: Pop,
                crossPlan: Union[str, np.ndarray] = "testcross",
                returnHybridPop: bool = False,
                simParam: Optional[SimParam] = None) -> Union[Pop, HybridPop]:
    """Convenience function with AlphaSimR-style parameter names"""
    return hybrid_cross(females, males, crossPlan, returnHybridPop, simParam)


def calcGCA(pop: Union[Pop, HybridPop], use: str = "pheno") -> Dict[str, np.ndarray]:
    """Convenience function with AlphaSimR-style parameter names"""
    return calc_gca(pop, use)


def setPhenoGCA(pop: Pop, testers: Pop, use: str = "pheno",
                h2: Optional[Union[float, List[float]]] = None,
                H2: Optional[Union[float, List[float]]] = None,
                varE: Optional[Union[float, List[float], np.ndarray]] = None,
                corE: Optional[np.ndarray] = None,
                reps: int = 1, fixEff: int = 1, p: Optional[float] = None,
                inbred: bool = False, onlyPheno: bool = False,
                simParam: Optional[SimParam] = None) -> Union[Pop, np.ndarray]:
    """Convenience function with AlphaSimR-style parameter names"""
    return set_pheno_gca(pop, testers, use, h2, H2, varE, corE, reps, fixEff, p, inbred, onlyPheno, simParam)


def self_pop(pop: Pop, n_progeny: int = 1,
             parents: Optional[List[int]] = None,
             keep_parents: bool = True,
             sim_param: Optional[SimParam] = None) -> Pop:
    """
    Create selfed progeny from each individual in a population.
    Mirrors AlphaSimR::self behavior for single-population objects.
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")

    if pop.ploidy % 2 != 0:
        raise ValueError("You can not self aneuploids")

    n_ind = pop.n_ind
    if parents is None:
        parent_idx = list(range(n_ind))
    else:
        parent_idx = [int(i) for i in parents]

    # Build cross plan: each selected individual selfed n_progeny times
    cross_plan = []
    for p in parent_idx:
        for _ in range(n_progeny):
            cross_plan.append((p, p))

    total_progeny = len(cross_plan)

    progeny_geno = []
    for chr_idx in range(pop.n_chr):
        n_packed = pop.geno[chr_idx].shape[0]
        chr_geno = np.zeros((n_packed, pop.ploidy, total_progeny), dtype=np.uint8)
        for prog_idx, (mother_idx, father_idx) in enumerate(cross_plan):
            for packed_idx in range(n_packed):
                for haplo in range(pop.ploidy):
                    # Randomly draw gametes from the two parental chromosome copies
                    if np.random.random() < 0.5:
                        chr_geno[packed_idx, haplo, prog_idx] = pop.geno[chr_idx][packed_idx, haplo, mother_idx]
                    else:
                        chr_geno[packed_idx, haplo, prog_idx] = pop.geno[chr_idx][packed_idx, haplo, father_idx]
        progeny_geno.append(chr_geno)

    # Parents for pedigree
    if keep_parents:
        mother = [pop.mother[p] for p, _ in cross_plan]
        father = [pop.father[p] for _, p in cross_plan]
    else:
        mother = [pop.id[p] for p, _ in cross_plan]
        father = [pop.id[p] for _, p in cross_plan]

    iid = list(range(sim_param._last_id + 1, sim_param._last_id + total_progeny + 1))
    ids = [str(i) for i in iid]

    progeny_pop = Pop(
        n_ind=total_progeny,
        n_chr=pop.n_chr,
        ploidy=pop.ploidy,
        n_loci=pop.n_loci,
        geno=progeny_geno,
        gen_map=pop.gen_map,
        centromere=pop.centromere,
        inbred=pop.inbred,
        id=ids,
        iid=iid,
        mother=mother,
        father=father,
        sex=["H"] * total_progeny,
        n_traits=pop.n_traits,
        gv=np.zeros((total_progeny, pop.n_traits)),
        pheno=np.zeros((total_progeny, pop.n_traits)),
        ebv=np.empty((total_progeny, 0)),
        gxe=[None] * pop.n_traits,
        fix_eff=[1] * total_progeny,
        misc={},
        misc_pop={}
    )

    if pop.n_traits > 0:
        progeny_pop.gv = _get_gv_index(progeny_pop, sim_param)
        progeny_pop.pheno = progeny_pop.gv.copy()

    sim_param._last_id += total_progeny

    return progeny_pop


def self(pop: Pop, nProgeny: int = 1,
         parents: Optional[List[int]] = None,
         keepParents: bool = True,
         simParam: Optional[SimParam] = None) -> Pop:
    """AlphaSimR-style self() convenience function"""
    return self_pop(pop, nProgeny, parents, keepParents, simParam)


def select_within_fam(pop: Pop, n_ind: int,
                      trait: Union[int, callable] = 1,
                      use: Union[str, callable] = "pheno",
                      sex: str = "B",
                      fam_type: str = "B",
                      select_top: bool = True,
                      return_pop: bool = True,
                      candidates: Optional[List[int]] = None,
                      sim_param: Optional[SimParam] = None,
                      **kwargs) -> Union[Pop, List[int]]:
    """
    Select a subset of individuals within each family.
    Mirrors AlphaSimR::selectWithinFam behavior for single Pop objects.
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")

    if n_ind < 0:
        raise ValueError("nInd must be >= 0")

    # Determine eligible individuals based on sex and optional candidates
    eligible = []
    for i in range(pop.n_ind):
        if sim_param._sexes == "no" or sex == "B" or pop.sex[i] == sex:
            eligible.append(i)

    if candidates is not None:
        eligible = [i for i in eligible if i in candidates]

    # Determine families
    fam_type = fam_type.upper()
    if fam_type == "B":
        families = [f"{pop.mother[i]}_{pop.father[i]}" for i in range(pop.n_ind)]
    elif fam_type == "F":
        families = pop.mother
    elif fam_type == "M":
        families = pop.father
    else:
        raise ValueError(f"famType={fam_type} is not a valid option")

    # Get selection criterion values
    if isinstance(use, str):
        if use == "pheno":
            values = pop.pheno[:, trait - 1] if isinstance(trait, int) else pop.pheno
        elif use == "gv":
            values = pop.gv[:, trait - 1] if isinstance(trait, int) else pop.gv
        elif use == "ebv":
            values = pop.ebv[:, trait - 1] if isinstance(trait, int) else pop.ebv
        elif use == "rand":
            values = np.random.random(pop.n_ind)
        else:
            raise ValueError(f"Unknown use criterion: {use}")
    else:
        values = use(pop, trait, **kwargs)

    if callable(trait):
        values = trait(values, **kwargs)

    if hasattr(values, "ndim") and values.ndim > 1:
        values = values.flatten()

    warn = False
    selected_indices: List[int] = []

    unique_fams = sorted(set(families))
    for fam in unique_fams:
        # indices belonging to this family
        fam_idx = [i for i, f in enumerate(families) if f == fam]
        fam_idx = [i for i in fam_idx if i in eligible]
        if not fam_idx:
            continue
        fam_vals = np.array([values[i] for i in fam_idx])
        order = np.argsort(fam_vals)
        if select_top:
            order = order[::-1]
        ordered_idx = [fam_idx[i] for i in order]
        if len(ordered_idx) < n_ind:
            warn = True
            take = ordered_idx
        else:
            take = ordered_idx[:n_ind]
        selected_indices.extend(take)

    if warn:
        # Mirror AlphaSimR behavior: warn but continue
        import warnings
        warnings.warn("One or more families are smaller than nInd")

    if return_pop:
        if not selected_indices:
            return Pop(
                n_ind=0,
                n_chr=pop.n_chr,
                ploidy=pop.ploidy,
                n_loci=pop.n_loci,
                geno=[],
                gen_map=pop.gen_map,
                centromere=pop.centromere,
                inbred=pop.inbred,
                id=[],
                iid=[],
                mother=[],
                father=[],
                sex=[],
                n_traits=pop.n_traits,
                gv=np.empty((0, pop.n_traits)),
                pheno=np.empty((0, pop.n_traits)),
                ebv=np.empty((0, 0)),
                gxe=[],
                fix_eff=[],
                misc={},
                misc_pop={}
            )

        return Pop(
            n_ind=len(selected_indices),
            n_chr=pop.n_chr,
            ploidy=pop.ploidy,
            n_loci=pop.n_loci,
            geno=[pop.geno[chr_idx][:, :, selected_indices] for chr_idx in range(pop.n_chr)],
            gen_map=pop.gen_map,
            centromere=pop.centromere,
            inbred=pop.inbred,
            id=[pop.id[i] for i in selected_indices],
            iid=[pop.iid[i] for i in selected_indices],
            mother=[pop.mother[i] for i in selected_indices],
            father=[pop.father[i] for i in selected_indices],
            sex=[pop.sex[i] for i in selected_indices],
            n_traits=pop.n_traits,
            gv=pop.gv[selected_indices, :],
            pheno=pop.pheno[selected_indices, :],
            ebv=pop.ebv[selected_indices, :],
            gxe=pop.gxe,
            fix_eff=[pop.fix_eff[i] for i in selected_indices],
            misc=pop.misc,
            misc_pop=pop.misc_pop
        )
    else:
        return selected_indices


def selectWithinFam(pop: Pop, nInd: int,
                    trait: Union[int, callable] = 1,
                    use: Union[str, callable] = "pheno",
                    sex: str = "B",
                    famType: str = "B",
                    selectTop: bool = True,
                    returnPop: bool = True,
                    candidates: Optional[List[int]] = None,
                    simParam: Optional[SimParam] = None,
                    **kwargs) -> Union[Pop, List[int]]:
    """AlphaSimR-style selectWithinFam() convenience function"""
    return select_within_fam(pop, nInd, trait, use, sex, famType,
                             selectTop, returnPop, candidates, simParam, **kwargs)


def make_dh(pop: Pop, n_dh: int = 1, use_female: bool = True,
            keep_parents: bool = True, sim_param: Optional[SimParam] = None) -> Pop:
    """
    Creates DH lines from each individual in a population.
    Only works with diploid individuals. For polyploids, use
    reduceGenome and doubleGenome.
    
    Parameters:
    -----------
    pop : Pop
        Population object
    n_dh : int, default=1
        Total number of DH lines per individual
    use_female : bool, default=True
        Should female recombination rates be used
    keep_parents : bool, default=True
        Should previous parents be used for mother and father
    sim_param : SimParam, optional
        Simulation parameters
    
    Returns:
    --------
    Pop
        New population of DH lines
    """
    if sim_param is None:
        raise ValueError("simParam must be provided")
    
    if pop.ploidy != 2:
        raise ValueError("Only works with diploids")
    
    n_ind = pop.n_ind
    total_dh = n_ind * n_dh
    
    # Select genetic map based on use_female
    if use_female:
        gen_map = sim_param._female_map
    else:
        if sim_param._male_map is None:
            gen_map = sim_param._female_map  # Fallback to female map
        else:
            gen_map = sim_param._male_map
    
    # Create DH genotypes
    dh_geno = []
    for chr_idx in range(pop.n_chr):
        n_packed = pop.geno[chr_idx].shape[0]
        n_loci_chr = pop.n_loci[chr_idx]
        chr_gen_map = gen_map[chr_idx]
        gen_len = chr_gen_map[-1] if len(chr_gen_map) > 0 else 0.0
        
        chr_dh_geno = np.zeros((n_packed, 2, total_dh), dtype=np.uint8)
        
        # For each individual, create nDH gametes and double them
        for ind_idx in range(n_ind):
            for dh_idx in range(n_dh):
                dh_global_idx = ind_idx * n_dh + dh_idx
                
                # Create a gamete through recombination simulation
                # Sample crossover positions based on genetic map (simplified Poisson process)
                # Using Kosambi mapping function parameter v from sim_param
                v = sim_param.v
                p = sim_param.p
                
                # Expected number of crossovers (simplified)
                # For Kosambi: mean = gen_len, variance depends on v
                # Simplified: sample crossovers as Poisson process
                if gen_len > 0:
                    # Sample crossover positions along genetic map
                    n_crossovers = np.random.poisson(gen_len)
                    if n_crossovers > 0:
                        crossover_positions = np.sort(np.random.uniform(0, gen_len, n_crossovers))
                        # Thin crossovers (keep 50% randomly)
                        crossover_positions = crossover_positions[np.random.random(n_crossovers) > 0.5]
                    else:
                        crossover_positions = np.array([])
                else:
                    crossover_positions = np.array([])
                
                # Initialize gamete with first chromosome
                gamete = pop.geno[chr_idx][:, 0, ind_idx].copy()
                
                # Apply crossovers: alternate between chromosomes at each crossover
                if len(crossover_positions) > 0:
                    # Start with chromosome 0, then alternate
                    current_chr = 0
                    prev_locus = 0
                    
                    for cross_pos in crossover_positions:
                        # Find which locus this crossover is at
                        cross_locus = np.searchsorted(chr_gen_map, cross_pos)
                        cross_locus = min(cross_locus, n_loci_chr - 1)
                        
                        # Copy from current chromosome up to this crossover
                        start_byte = prev_locus // 8
                        end_byte = cross_locus // 8
                        
                        for byte_idx in range(start_byte, min(end_byte + 1, n_packed)):
                            gamete[byte_idx] = pop.geno[chr_idx][byte_idx, current_chr, ind_idx]
                        
                        # Switch to other chromosome for next segment
                        current_chr = 1 - current_chr
                        prev_locus = cross_locus
                    
                    # Copy remaining from current chromosome
                    start_byte = prev_locus // 8
                    for byte_idx in range(start_byte, n_packed):
                        gamete[byte_idx] = pop.geno[chr_idx][byte_idx, current_chr, ind_idx]
                
                # Double the gamete to create diploid DH (both copies identical)
                chr_dh_geno[:, 0, dh_global_idx] = gamete
                chr_dh_geno[:, 1, dh_global_idx] = gamete
        
        dh_geno.append(chr_dh_geno)
    
    # Create MapPop for DH lines
    raw_pop = MapPop(
        n_ind=total_dh,
        n_chr=pop.n_chr,
        ploidy=pop.ploidy,
        n_loci=pop.n_loci,
        geno=dh_geno,
        gen_map=pop.gen_map,
        centromere=pop.centromere,
        inbred=True  # DH lines are inbred
    )
    
    # Set parent information
    if keep_parents:
        mother = [pop.mother[i] for i in range(n_ind) for _ in range(n_dh)]
        father = [pop.father[i] for i in range(n_ind) for _ in range(n_dh)]
        i_mother = [pop.iid[i] for i in range(n_ind) for _ in range(n_dh)]
        i_father = [pop.iid[i] for i in range(n_ind) for _ in range(n_dh)]
    else:
        mother = [pop.id[i] for i in range(n_ind) for _ in range(n_dh)]
        father = [pop.id[i] for i in range(n_ind) for _ in range(n_dh)]
        i_mother = [pop.iid[i] for i in range(n_ind) for _ in range(n_dh)]
        i_father = [pop.iid[i] for i in range(n_ind) for _ in range(n_dh)]
    
    # Generate IDs for DH lines
    iid = list(range(sim_param._last_id + 1, sim_param._last_id + total_dh + 1))
    id_list = [str(i) for i in iid]
    
    # Create Pop object
    dh_pop = Pop(
        n_ind=total_dh,
        n_chr=pop.n_chr,
        ploidy=pop.ploidy,
        n_loci=pop.n_loci,
        geno=dh_geno,
        gen_map=pop.gen_map,
        centromere=pop.centromere,
        inbred=True,
        id=id_list,
        iid=iid,
        mother=mother,
        father=father,
        sex=["H"] * total_dh,
        n_traits=pop.n_traits,
        gv=np.zeros((total_dh, pop.n_traits)),
        pheno=np.zeros((total_dh, pop.n_traits)),
        ebv=np.empty((total_dh, 0)),
        gxe=[None] * pop.n_traits,
        fix_eff=[1] * total_dh,
        misc={},
        misc_pop={}
    )
    
    # Calculate genetic values for DH lines
    if pop.n_traits > 0:
        dh_pop.gv = _get_gv_index(raw_pop, sim_param)
        dh_pop.pheno = dh_pop.gv.copy()
    
    # Update last ID
    sim_param._last_id += total_dh
    
    return dh_pop


def makeDH(pop: Pop, nDH: int = 1, useFemale: bool = True,
           keepParents: bool = True, simParam: Optional[SimParam] = None) -> Pop:
    """AlphaSimR-style makeDH() convenience function"""
    return make_dh(pop, nDH, useFemale, keepParents, simParam)