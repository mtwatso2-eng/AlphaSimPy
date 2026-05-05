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

   **Check if RcppArmadillo is already installed:**
   ```bash
   R --slave -e "if (require('RcppArmadillo', quietly=TRUE)) cat('RcppArmadillo is installed\n') else cat('RcppArmadillo is NOT installed\n')"
   ```

   **If not installed, use one of these methods:**
   
   Option A (using Rscript):
   ```bash
   Rscript -e "install.packages('RcppArmadillo', repos='https://cran.r-project.org')"
   ```
   
   Option B (using full path if Rscript not in PATH):
   ```bash
   /usr/local/bin/Rscript -e "install.packages('RcppArmadillo', repos='https://cran.r-project.org')"
   ```
   
   Option C (using R directly):
   ```bash
   R --slave -e "install.packages('RcppArmadillo', repos='https://cran.r-project.org')"
   ```
   
   **Note:** The build script will automatically find RcppArmadillo if it's installed, so you can skip this step if it's already available.

3. Build the extension:
```bash
python setup_alphasimpy.py build_ext --inplace
```

Or install in development mode (if using a pyproject.toml that references this):
```bash
pip install -e .
```

**Note:** Use `setup_alphasimpy.py` for building the C++ bindings. The main `setup.py` in this project is used for other purposes.

## Troubleshooting

### R/RcppArmadillo not found
If the build fails to find R or RcppArmadillo:
- Ensure R is in your PATH (or use full path: `/usr/local/bin/R`)
- Check if RcppArmadillo is installed: `R --slave -e "require('RcppArmadillo')"`
- If Rscript doesn't work, use `R --slave -e` instead
- Set `R_HOME` environment variable if needed: `export R_HOME=$(R RHOME)`
- Install R development packages (e.g., `r-base-dev` on Ubuntu)
- The setup script automatically searches for RcppArmadillo via R, so R must be accessible

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
