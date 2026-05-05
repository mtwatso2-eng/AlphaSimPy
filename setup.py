#!/usr/bin/env python3
"""Setuptools build for AlphaSimPy package and C++ extension."""

import os
import sys
from pathlib import Path

from setuptools import Extension, setup
from setuptools.command.build_ext import build_ext

try:
    from pybind11.setup_helpers import Pybind11Extension
    HAS_PYBIND11_HELPER = True
except ImportError:
    HAS_PYBIND11_HELPER = False


def get_explicit_alphasimr_src():
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
            "Missing AlphaSimR C++ path. Provide --alphasimr-src=/path/to/src or set ALPHASIMR_SRC."
        )

    src_path = Path(raw_path).expanduser().resolve()
    if not src_path.exists() or not src_path.is_dir():
        raise RuntimeError(f"Invalid AlphaSimR source directory: {src_path}")
    return src_path


def get_boost_include():
    env_path = os.environ.get("BOOST_INCLUDE_DIR")
    if env_path and os.path.exists(os.path.join(env_path, "boost")):
        return env_path
    for path in ["/usr/local/include", "/opt/homebrew/include", "/opt/local/include", "/usr/include"]:
        if os.path.exists(os.path.join(path, "boost")):
            return path
    return None


def get_armadillo_include():
    env_path = os.environ.get("ARMADILLO_INCLUDE_DIR")
    if env_path and os.path.exists(os.path.join(env_path, "armadillo")):
        return env_path
    for path in ["/usr/local/include", "/opt/homebrew/include", "/opt/local/include", "/usr/include"]:
        if os.path.exists(os.path.join(path, "armadillo")):
            return path
    return None


class BuildExt(build_ext):
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
                "Boost headers not found. Install system Boost or set BOOST_INCLUDE_DIR."
            )

        super().build_extension(ext)


ALPHASIMR_SRC = get_explicit_alphasimr_src()
ALPHASIMPY_SRC = Path(__file__).parent / "src"

alphasimr_sources = []
for f in ["algorithm.cpp", "datastructures.cpp", "simulator.cpp", "misc_standalone.cpp"]:
    p = ALPHASIMR_SRC / f
    if p.exists():
        alphasimr_sources.append(str(p))

include_dirs = [str(ALPHASIMR_SRC), str(ALPHASIMPY_SRC)]
boost_inc = get_boost_include()
if boost_inc:
    include_dirs.append(boost_inc)
armadillo_inc = get_armadillo_include()
if armadillo_inc:
    include_dirs.append(armadillo_inc)

try:
    import pybind11
    include_dirs.insert(0, pybind11.get_include())
except ImportError:
    pass

ext_name = "alphasimpy._alphasimpy_cpp"
ext_sources = [str(ALPHASIMPY_SRC / "alphasimpy_bindings.cpp")] + alphasimr_sources

if HAS_PYBIND11_HELPER:
    ext_modules = [
        Pybind11Extension(
            ext_name,
            sources=ext_sources,
            include_dirs=include_dirs,
            cxx_std=14,
            extra_compile_args=["-DARMA_64BIT_WORD", "-DBOOST_DISABLE_ASSERTS"]
            + (["-fopenmp"] if sys.platform != "darwin" else []),
            extra_link_args=(["-fopenmp"] if sys.platform != "darwin" else ["-undefined", "dynamic_lookup"]),
        )
    ]
else:
    ext_modules = [
        Extension(
            ext_name,
            sources=ext_sources,
            include_dirs=include_dirs,
            extra_compile_args=["-std=c++14", "-DARMA_64BIT_WORD", "-DBOOST_DISABLE_ASSERTS"]
            + (["-fopenmp"] if sys.platform != "darwin" else []),
            extra_link_args=(["-fopenmp"] if sys.platform != "darwin" else ["-undefined", "dynamic_lookup"]),
        )
    ]


setup(
    name="alphasimpy",
    version="0.3.0",
    description="AlphaSimPy: breeding simulation with required C++/pybind11 extension",
    packages=["alphasimpy"],
    package_dir={"": "src"},
    py_modules=["AlphaSimPy"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": BuildExt},
    install_requires=["numpy>=1.19.0", "pybind11>=2.6.0"],
    python_requires=">=3.9",
)