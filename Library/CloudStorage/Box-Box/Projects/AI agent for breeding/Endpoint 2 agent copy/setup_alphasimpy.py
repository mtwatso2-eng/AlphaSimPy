#!/usr/bin/env python3
"""
Setup script for building AlphaSimPy C++ bindings.

Builds the _alphasimpy_cpp extension module that wraps AlphaSimR's MaCS
simulation code. Uses a standalone C++ toolchain (no R dependency).

Usage:
    python setup_alphasimpy.py build_ext --inplace
    # or
    pip install -e .  (if using pyproject.toml that points to this)
"""

import os
import sys
from pathlib import Path

from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext

# Try to use pybind11's setup helpers
try:
    from pybind11.setup_helpers import Pybind11Extension, build_ext as PybindBuildExt
    from pybind11 import get_cmake_dir
    HAS_PYBIND11_HELPER = True
except ImportError:
    HAS_PYBIND11_HELPER = False


def get_explicit_alphasimr_src():
    """
    Resolve AlphaSimR C++ source directory.

    Resolution order:
    - Bundled folder: ./src/alphasimr_cpp
    - CLI: --alphasimr-src=/path/to/AlphaSimR-master/src
    - ENV: ALPHASIMR_SRC=/path/to/AlphaSimR-master/src
    """
    bundled = Path(__file__).parent / "src" / "alphasimr_cpp"
    if bundled.exists() and bundled.is_dir():
        raw_path = str(bundled)
    else:
        raw_path = None

    cli_value = None
    for arg in list(sys.argv):
        if arg.startswith("--alphasimr-src="):
            cli_value = arg.split("=", 1)[1]
            sys.argv.remove(arg)
            break

    raw_path = cli_value or os.environ.get("ALPHASIMR_SRC") or raw_path
    if not raw_path:
        raise RuntimeError(
            "Missing AlphaSimR C++ path. Provide --alphasimr-src=/path/to/AlphaSimR-master/src, "
            "set ALPHASIMR_SRC, or include bundled sources at ./src/alphasimr_cpp."
        )

    src_path = Path(raw_path).expanduser().resolve()
    if not src_path.exists() or not src_path.is_dir():
        raise RuntimeError(f"Invalid AlphaSimR source directory: {src_path}")

    required_files = ["algorithm.cpp", "datastructures.cpp", "simulator.cpp", "alphasimr.h"]
    missing = [name for name in required_files if not (src_path / name).exists()]
    if missing:
        raise RuntimeError(
            f"AlphaSimR source directory is missing required files: {missing}. "
            f"Expected folder like .../AlphaSimR-master/src, got: {src_path}"
        )

    return src_path


def get_boost_include():
    """Find Boost include directory from system package managers."""
    env_path = os.environ.get("BOOST_INCLUDE_DIR")
    if env_path and os.path.exists(os.path.join(env_path, "boost")):
        return env_path
    for path in [
        "/usr/local/include",
        "/opt/homebrew/include",
        "/opt/local/include",
        "/usr/include",
    ]:
        if os.path.exists(os.path.join(path, "boost")):
            return path
    return None


def get_armadillo_include():
    """Find Armadillo headers in common system locations."""
    env_path = os.environ.get("ARMADILLO_INCLUDE_DIR")
    if env_path and os.path.exists(os.path.join(env_path, "armadillo")):
        return env_path
    for path in [
        "/usr/local/include",
        "/opt/homebrew/include",
        "/opt/local/include",
        "/usr/include",
    ]:
        if os.path.exists(os.path.join(path, "armadillo")):
            return path
    return None


class BuildExt(build_ext):
    """Custom build_ext for standalone (no-R) simulator build."""

    def build_extension(self, ext):
        armadillo_inc = get_armadillo_include()
        if armadillo_inc and armadillo_inc not in ext.include_dirs:
            ext.include_dirs.append(armadillo_inc)
        if "armadillo" not in ext.libraries:
            ext.libraries.append("armadillo")

        boost_inc = get_boost_include()
        if boost_inc and boost_inc not in ext.include_dirs:
            ext.include_dirs.append(boost_inc)
        elif not boost_inc:
            raise RuntimeError(
                "Boost headers not found. Install via: "
                "brew install boost (macOS), apt install libboost-all-dev (Linux), "
                "or set BOOST_INCLUDE_DIR."
            )

        # Ensure we are not accidentally building against R-provided headers.
        for inc in ext.include_dirs:
            norm = os.path.abspath(inc)
            if "/R/" in norm or norm.endswith("/R"):
                raise RuntimeError(
                    f"R include path detected in no-R build: {norm}. "
                    "Install standalone Boost/Armadillo and remove R include paths."
                )

        super().build_extension(ext)


# AlphaSimR source files needed for MaCS (excluding RcppExports)
ALPHASIMR_SRC = get_explicit_alphasimr_src()
ALPHASIMPY_SRC = Path(__file__).parent / "src"

alphasimr_sources = []
if ALPHASIMR_SRC.exists():
    for f in ["algorithm.cpp", "datastructures.cpp"]:
        p = ALPHASIMR_SRC / f
        if p.exists():
            alphasimr_sources.append(str(p))
    # simulator.cpp will be patched at build time
    sim_path = ALPHASIMR_SRC / "simulator.cpp"
    if sim_path.exists():
        alphasimr_sources.append(str(sim_path))
    standalone_misc = ALPHASIMR_SRC / "misc_standalone.cpp"
    if standalone_misc.exists():
        alphasimr_sources.append(str(standalone_misc))

include_dirs = [
    str(ALPHASIMR_SRC) if ALPHASIMR_SRC.exists() else ".",
    str(ALPHASIMPY_SRC),
]

boost_inc = get_boost_include()
if boost_inc:
    include_dirs.append(boost_inc)
else:
    print(
        "Warning: Boost headers not found. "
        "Set BOOST_INCLUDE_DIR or install Boost (e.g. brew install boost)."
    )

armadillo_inc = get_armadillo_include()
if armadillo_inc:
    include_dirs.append(armadillo_inc)
else:
    print(
        "Warning: Armadillo headers not found. "
        "Set ARMADILLO_INCLUDE_DIR or install Armadillo (e.g. brew install armadillo)."
    )

# Add pybind11 include for non-helper case
try:
    import pybind11
    include_dirs.insert(0, pybind11.get_include())
except ImportError:
    pass

ext_modules = []
if HAS_PYBIND11_HELPER:
    ext_modules = [
        Pybind11Extension(
            "alphasimpy._alphasimpy_cpp",
            sources=[str(ALPHASIMPY_SRC / "alphasimpy_bindings.cpp")] + alphasimr_sources,
            include_dirs=include_dirs,
            cxx_std=14,
            extra_compile_args=[
                "-DARMA_64BIT_WORD",
                "-DBOOST_DISABLE_ASSERTS",
            ]
            + (["-fopenmp"] if sys.platform != "darwin" else []),
            extra_link_args=(
                ["-fopenmp"]
                if sys.platform != "darwin"
                else ["-undefined", "dynamic_lookup"]
            ),
        )
    ]
else:
    ext_modules = [
        Extension(
            "alphasimpy._alphasimpy_cpp",
            sources=[str(ALPHASIMPY_SRC / "alphasimpy_bindings.cpp")] + alphasimr_sources,
            include_dirs=include_dirs,
            extra_compile_args=[
                "-std=c++14",
                "-DARMA_64BIT_WORD",
                "-DBOOST_DISABLE_ASSERTS",
            ]
            + (["-fopenmp"] if sys.platform != "darwin" else []),
            extra_link_args=(
                ["-fopenmp"] if sys.platform != "darwin" else ["-undefined", "dynamic_lookup"]
            ),
        )
    ]

setup(
    name="alphasimpy",
    version="0.1.8",
    description="AlphaSimPy with C++ bindings for MaCS simulation",
    packages=["alphasimpy"],
    package_dir={"": "src"},
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
    install_requires=[
        "numpy>=1.19.0",
        "pybind11>=2.6.0",
    ],
)
