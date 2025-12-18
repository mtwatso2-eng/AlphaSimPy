# Building AlphaSimPy with C++ Bindings

AlphaSimPy now uses PyBind11 to call the original C++ code from AlphaSimR instead of Python translations.

## Prerequisites

1. **Python 3.7+** with development headers
2. **C++ compiler** with C++11 support (g++, clang++, etc.)
3. **R and RcppArmadillo** - The C++ code uses RcppArmadillo which requires:
   - R installed and accessible
   - RcppArmadillo package installed in R
   - R development headers
4. **Boost libraries** (for string splitting)
5. **Armadillo** (included via RcppArmadillo)

## Building the Extension

1. Install Python dependencies:
```bash
pip install -r requirements.txt
```

2. Ensure R is installed and RcppArmadillo is available:
```bash
Rscript -e "install.packages('RcppArmadillo')"
```

3. Build the extension:
```bash
python setup.py build_ext --inplace
```

Or install in development mode:
```bash
pip install -e .
```

## Troubleshooting

### R/RcppArmadillo not found
If the build fails to find R or RcppArmadillo:
- Ensure R is in your PATH
- Set `R_HOME` environment variable if needed
- Install R development packages (e.g., `r-base-dev` on Ubuntu)

### Boost not found
Install Boost development libraries:
- Ubuntu/Debian: `sudo apt-get install libboost-dev`
- macOS: `brew install boost`
- Or compile from source

### Compilation errors
- Ensure all C++ source files are in the `src/` directory
- Check that all header files are accessible
- Verify C++11 support is enabled

## Usage

After building, the Python code will automatically use the C++ bindings:

```python
import AlphaSimPy
founder_pop = AlphaSimPy.runMacs(nInd=100, nChr=1)
```

If the C++ bindings are not available, the code will raise an error (the Python fallback has been removed as requested).
