# AlphaSimPy

Python breeding simulation aligned with [AlphaSimR](https://github.com/gaynorr/AlphaSimR). The package builds an optional **C++ extension** (pybind11 + AlphaSimR-derived sources under `src/alphasimr_cpp/`) for faster coalescent-style founder generation and related paths. If the extension is missing, many routines **fall back to pure Python** (you will see a one-time import warning).

## Requirements

| | |
|---|---|
| **Python** | 3.9+ (see `pyproject.toml`) |
| **Runtime** | `numpy`, `pybind11` |
| **Build (extension)** | C++14-capable toolchain, **Boost headers**, **Armadillo** (headers + linked library) |

Bundled C++ sources live in `src/alphasimr_cpp/`. You do **not** need a separate AlphaSimR checkout unless you override the source path (see below).

## Quick install (from a clone)

### 1. System libraries

**macOS** (recommended):

```bash
xcode-select --install   # if you do not already have Command Line Tools
brew install boost armadillo
```

**Ubuntu / Debian**:

```bash
sudo apt-get update
sudo apt-get install -y build-essential \
  libboost-all-dev libarmadillo-dev
```

**Fedora-style**:

```bash
sudo dnf install gcc-c++ boost-devel armadillo-devel
```

If headers are installed in a non-standard location:

```bash
export BOOST_INCLUDE_DIR=/path/to/boost/parent       # directory containing boost/
export ARMADILLO_INCLUDE_DIR=/path/to/include         # directory containing armadillo/
```

### 2. Python environment and package

From the repository root:

```bash
python -m venv .venv
source .venv/bin/activate   # Windows: .venv\Scripts\activate
pip install -U pip setuptools wheel

# editable install while developing (builds extension in-place)
pip install -e .

# optional: notebooks, plotting, etc.
pip install -r requirements.txt
```

To force a fresh build of the extension after C++ changes:

```bash
pip install -e . --no-build-isolation --force-reinstall
```

### Advanced: alternate AlphaSimR C++ tree

Only if you are **not** using the bundled `src/alphasimr_cpp`, point `setup.py` at another AlphaSimR C++ directory:

```bash
export ALPHASIMR_SRC=/absolute/path/to/AlphaSimR/src
pip install -e .
```

If you invoke `setup.py` directly, you can instead pass **`--alphasimr-src=/path`**.

Tagged releases may ship prebuilt wheels from `.github/workflows/wheels.yml` (Linux/macOS).

## Smoke test

```bash
pytest -q
```

(`pyproject.toml` sets `pythonpath = ["src"]` for pytest when run from the repo root.)

```python
from alphasimpy import MapPop, SimParam, new_pop, run_macs

founder_pop = run_macs(n_ind=10, n_chr=1, seg_sites=100)
SP = SimParam(founder_pop)
SP.addTraitA(5)
SP.addSnpChip(10)
pop = new_pop(founder_pop, sim_param=SP)
print(pop.n_ind, pop.n_traits)
```

You can also `import AlphaSimPy` — the canonical package namespace **`alphasimpy`** re-exports the same symbols.

## Usage overview

### `run_macs`

Mirrors AlphaSimR `runMacs` semantics with Python-style names:

- **`n_ind`**, **`n_chr`**, **`seg_sites`**, **`inbred`**, **`species`** (`GENERIC`, `CATTLE`, `WHEAT`, `MAIZE`), **`split`**, **`ploidy`**, **`manual_command`** / **`manual_gen_len`**, **`n_threads`**

### Species models

- **GENERIC** — default population history  
- **CATTLE** — Macleod et al.-style demographic history  
- **WHEAT**, **MAIZE** — crop-specific presets  

### `SimParam` examples

```python
SP = SimParam(founder_pop)
SP.addTraitA(nQtlPerChr=5, mean=0, var=1)
SP.addTraitAD(nQtlPerChr=5, meanDD=0.5)
SP.addSnpChip(nSnpPerChr=100)
SP.setTrackPed(True)
SP.setTrackRec(True)
SP.setSexes("yes_sys")
```

### `Pop`

Build from founders with **`new_pop`**, **`new_empty_pop`**, **`new_multi_pop`**, etc. Main attributes include **`id`**, **`gv`**, **`pheno`**, **`ebv`**, pedigree fields, **`misc`** / **`misc_pop`**.

See `src/AlphaSimPy.py` and **`alphasimpy.__all__`** for the full public API (`run_macs`, crossing, selection, `merge_pops`, `edit_genome`, …).

### `MapPop` (datasummary)

Population with genetic map metadata: chromosome counts, packed **`geno`**, **`gen_map`**, **`centromere`**, **ploidy**, **inbred** flag.

## Troubleshooting builds

| Symptom | What to try |
|---------|-------------|
| `Boost headers not found` | Install Boost dev packages or set **`BOOST_INCLUDE_DIR`**. |
| Armadillo not found | Install `libarmadillo-dev` (Debian) / `brew install armadillo` / set **`ARMADILLO_INCLUDE_DIR`**. |
| Link errors for `armadillo` | Ensure the Armadillo **library** is installed, not headers only. |
| `WARNING: C++ bindings not available` | Extension failed to load; install rebuilt with steps above or use Python fallback (slower). |
| Wrong Python | Use `python -m pip install -e .` from the **same** interpreter you run. |

On **macOS**, OpenMP is disabled in `setup.py` for compatibility; Linux builds may use `-fopenmp` where supported.

## License

This project draws on ideas and C++ lineage from AlphaSimR; see AlphaSimR’s license for the original R package. Add or follow a `LICENSE` file in this repository for the Python port itself, if you redistribute it.
