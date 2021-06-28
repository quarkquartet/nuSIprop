# In this file, we declare the C++ objects and methods we want to import
cdef extern from "nuSIprop.hpp" namespace "nuSIprop":
    cdef cppclass calculate_flux:
        calculate_flux() except +
        calculate_flux(double, double, double, double, 
		       double,
		       bint, bint, bint,
		       int, double, double,
		       double, int, bint) except + # If exceptions are raised at construction, convert them to python exceptions
        # Methods and variables to import
        void evolve()
        double get_flux(int, int)
        double get_flux_fla(int, int)
        double get_energy(int)
        int get_N_bins_E()
        double check_energy_conservation()

        double mphi, g, mntot, si, norm
