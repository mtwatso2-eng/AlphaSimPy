# AlphaSimPy

Python clone of **[AlphaSimR](https://github.com/gaynorr/AlphaSimR)** for genomic breeding simulation. The package **requires** a compiled C extension (pybind11, Boost, Armadillo); `pip install …` succeeds only once that extension is available via a matching wheel or a local build.

## Installation

Use published wheels when they match your OS and Python (no compiler needed):

```bash
pip install alphasimpy
```

Otherwise clone the repo and build from source (**Python ≥3.9**). Install system dependencies, then editable-install from the repo root:

**macOS**

```bash
xcode-select --install  # Command Line Tools, if missing
brew install boost armadillo
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

**Ubuntu / Debian**

```bash
sudo apt-get install -y build-essential libboost-all-dev libarmadillo-dev
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

**Fedora**

```bash
sudo dnf install gcc-c++ boost-devel armadillo-devel
python -m venv .venv && source .venv/bin/activate
pip install -U pip setuptools wheel
pip install -e .
```

Non-default header locations:

```bash
export BOOST_INCLUDE_DIR=/path/to/parent/of/boost
export ARMADILLO_INCLUDE_DIR=/path/to/parent/of/armadillo
```

Extras (notebooks, etc.):

```bash
pip install -r requirements.txt
```

Rebuild the extension after C++ changes:

```bash
pip install -e . --no-build-isolation --force-reinstall
```
