import os
from distutils.core import Extension, setup

import numpy as np
from Cython.Build import cythonize

os.environ["CC"] = "/opt/homebrew/Cellar/gcc/14.1.0_1/bin/gcc-14"
os.environ["CXX"] = "/opt/homebrew/Cellar/gcc/14.1.0_1/bin/g++-14"

gsl_include_dir = "/opt/homebrew/opt/gsl/include"
gsl_lib_dir = "/opt/homebrew/opt/gsl/lib"

polylogarithm_sources = [
    "/Users/quarkquartet/Dropbox/Research_Project/2024-1-Massless_Neutrino_DSNB/02-Analysis/nuSIprop/polylogarithm/src/cpp/Li2.cpp",
    "/Users/quarkquartet/Dropbox/Research_Project/2024-1-Massless_Neutrino_DSNB/02-Analysis/nuSIprop/polylogarithm/src/cpp/Li3.cpp",  # Add more source files as needed
]

sources = ["nuSIprop.pyx"] + polylogarithm_sources

extensions = [
    Extension(
        "nuSIprop",
        sources=sources,
        include_dirs=[np.get_include(), gsl_include_dir],
        libraries=["gsl", "gslcblas", "m"],
        library_dirs=[gsl_lib_dir],
        extra_compile_args=["-O3", "-std=gnu++11", "-v"],  # Enable verbose output
        extra_link_args=["-lgsl", "-lgslcblas", "-lm", "-v"],  # Enable verbose output
    )
]

setup(
    ext_modules=cythonize(extensions),
)
