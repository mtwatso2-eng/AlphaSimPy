from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup, Extension
import pybind11
import os

# Get the source directory
src_dir = "src"

# Collect only the C++ source files needed for MaCS simulation
# We need simulator.cpp, algorithm.cpp, misc.cpp, and datastructures.cpp
cpp_sources = [
    os.path.join(src_dir, "alphasimpy_bindings.cpp"),
    os.path.join(src_dir, "simulator.cpp"),
    os.path.join(src_dir, "algorithm.cpp"),
    os.path.join(src_dir, "misc.cpp"),
    os.path.join(src_dir, "datastructures.cpp"),
]

# Try to find Armadillo
import subprocess
armadillo_cflags = ""
armadillo_libs = ""
try:
    result = subprocess.run(['pkg-config', '--cflags', 'armadillo'], 
                          capture_output=True, text=True)
    if result.returncode == 0:
        armadillo_cflags = result.stdout.strip()
    result = subprocess.run(['pkg-config', '--libs', 'armadillo'],
                          capture_output=True, text=True)
    if result.returncode == 0:
        armadillo_libs = result.stdout.strip()
except:
    pass

# Include directories
include_dirs = [
    src_dir,
    pybind11.get_include(),
]

# Add Armadillo include flags if found
if armadillo_cflags:
    import shlex
    cflags_list = shlex.split(armadillo_cflags)
    for flag in cflags_list:
        if flag.startswith('-I'):
            include_dirs.append(flag[2:])


import os
import subprocess
import shlex

# Compiler flags
compile_args = [
    "-std=c++14",  # Armadillo requires C++14
    "-O3",
    "-DARMA_DONT_USE_WRAPPER",
    "-DALPHASIMPY_PYTHON",
]

# Add Armadillo compile flags if found via pkg-config (optional)
if armadillo_cflags:
    cflags_list = shlex.split(armadillo_cflags)
    for flag in cflags_list:
        if flag.startswith('-I'):
            include_dirs.append(flag[2:])
        else:
            compile_args.append(flag)

# Linker flags
link_args = []

# 1) Explicitly add the conda env lib dir
conda_prefix = os.environ.get("CONDA_PREFIX")
if conda_prefix:
    link_args.append(f"-L{os.path.join(conda_prefix, 'lib')}")

# 2) Add Armadillo library flag
if armadillo_libs:
    libs_list = shlex.split(armadillo_libs)
    link_args.extend(libs_list)
else:
    link_args.append("-larmadillo")

# Define the extension module
ext_modules = [
    Pybind11Extension(
        "_alphasimpy_cpp",
        cpp_sources,
        include_dirs=include_dirs,
        language="c++",
        cxx_std=14,
        extra_compile_args=compile_args,
        extra_link_args=link_args,
    ),
]

setup(
    name="alphasimpy",
    version="0.1.0",
    author="AlphaSimPy",
    description="AlphaSimPy with PyBind11 C++ bindings",
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    zip_safe=False,
    python_requires=">=3.7",
    install_requires=[
        "pybind11>=2.6.0",
        "numpy>=1.19.0",
    ],
)
