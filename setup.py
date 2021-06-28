from distutils.core import setup, Extension
from Cython.Build import cythonize

extensions = [Extension("nuSIprop", ["nuSIprop.pyx"],
                        #include_dirs=[...],
                        libraries=["gsl", "gslcblas"] ,
                        # library_dirs=[...],
                        extra_compile_args=["-O3", "-std=gnu++11"] 
                        )]
setup(
    ext_modules = cythonize(extensions),
)
