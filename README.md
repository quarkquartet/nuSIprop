# nuSIprop
This code simulates propagation of astrophysical neutrinos under the presence of neutrino self-interactions. Neutrinos are assumed to follow a power law in energy, and their sources are assumed to be distributed in redshift according to the star formation rate. We also assume that all mass eigenstates are equally produced at sources, and that interactions are diagonal in flavor space. For the neutrino mixing parameters, we assume the [NuFIT 5.0](http://www.nu-fit.org/?q=node/228) best fit values. All these assumptions can be lifted through simple modifications of the code.

For information on the physics, see the companion paper. Please cite it if you use this code. For details on the numerical algorithm, check the file `Details.pdf`.

## Prerequisites
To run nuSIprop, you will need
* `gcc`, or other C++ compiler supporting the GNU C++ extensions.
* `GSL`

Additionally, the Python interface requires
* `python3`
* `Cython`
* `numpy`
* `scipy`

## Basic usage
We provide both C++ and Python interfaces. The file `test.cpp` contains a minimal working example of C++ usage. Compile it as

```g++ -std=gnu++11 -O3 -o test.out test.cpp -lgsl -lgslcblas```

After running it, execute

```./test.out```

This will output the neutrino flux for each flavor as a function of energy assuming self-interactions.

Alternatively, one can use the Python interface. First compile the Python library by running

```python3 setup.py build_ext --inplace```

and then run the example with

```python3 test.py```

Both test files contain, at the end, additional usage tips and a detailed description of the available functions.

## Advanced usage: double scalar production
The (double) integrals associated to scattering through on-shell production of two mediators cannot be done analytically. 
Because of this, we have chosen to pre-compute them, store them in tables, and then interpolate these tables. This channel is numerically subleading and turned off by default, but it can be activated using the `phiphi` option in the object constructor.

To generate the interpolating tables, the required code is in the folder `xsec`. One has to first compile `funcs.c` by running

```gcc -shared -o funcs.so -fPIC funcs.c```

And then the tables can be generated by running the python script

```python3 tables_phiphi.py```

For improved access speed and file size, the `text_to_binary.cpp` script converts the `alpha_phiphi.dat` and `alphatilde_phiphi.dat` text files into binary files. Any C++ compiler should be able to compile it.

This procedure should generate two files, `alpha_phiphi.bin` and `alphatilde_phiphi.bin`, that the main code uses to include double scalar production. These two `.bin` files are also available upon request.
