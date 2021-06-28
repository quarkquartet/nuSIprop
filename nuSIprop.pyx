# In this file, we create a python-style class wrapping the C++ code
# Just run python setup.py build_ext --inplace

# IMPORTANT: tell the compiler that this is c++
# distutils: language = c++ 

from nuSIprop cimport calculate_flux
import numpy as np
import warnings
import scipy.interpolate as interp

cdef class pyprop:
    """Class that evolves an astrophysical neutrino flux assuming self-interactions

        Mandatory parameters:
        mphi   ---- Mediator mass [eV]
        g      ---- Coupling. The interaction Lagrangian is
                    - (1/2) g  \bar{\psi}\psi \phi,
                    with \psi and \phi the neutrino and scalar fields,
                    respectively; and the 1/2 factor is for Majorana
                    neutrinos.
        mntot  ---- Sum of neutrino masses [eV]
        si     ---- Spectral index of the injected flux

        The energy dependence of the injected flux is given by
        f(E) = (E/100 TeV)^{-si }

        Optional parameters:

        norm            ---- Normalization of the injected flux. Default: 1
                             This will be the free-streaming flavor-summed neutrino+antineutrino flux arriving at Earth at E=100 TeV
        majorana        ---- Whether to consider Majorana (true) or Dirac (false) neutrinos. Default: true
        non_resonant    ---- Whether to consider non s- scattering channels (true) or not (false). Default: true
        normal_ordering ---- Whether to consider normal mass ordering (true) or inverted (false). Default: true
        N_bins_E        ---- Number of energy bins. Default: 300
        lEmin           ---- log_10(smallest energy considered / eV). Default: 12
        lEmax           ---- log_10(largest energy considered / eV). Default: 17
        zmax            ---- Largest redsfhit. Default: 5
        flav            ---- Flavor of interacting neutrinos: 0=e, 1=mu, 2=tau. Default: 2
        phiphi          ---- Whether to include double scalar production (true) or not (false).
                             This requires having the cross-section tables. Default: true
    """
    cdef calculate_flux my_object  # Hold a C++ instance which we're wrapping

    cdef bint evolved # Wheter the evolve() function has been called

    def __cinit__(self,
                  double mphi, double g, double mntot, double si,
		  double norm = 1,
		  bint majorana = True, bint non_resonant = True, bint normal_ordering = True,
		  int N_bins_E = 300, double lEmin = 12.0, double lEmax = 17.0,
		  double zmax = 5.0, int flav=2, bint phiphi = True):
        self.my_object = calculate_flux(mphi, g, mntot, si,
                                        norm,
                                        majorana, non_resonant, normal_ordering,
                                        N_bins_E, lEmin, lEmax,
                                        zmax, flav, phiphi)
        self.evolved = False

    def set_parameters(self, mphi = None, g = None, mntot = None, si = None, norm = None):
        """Modify the physics parameters
        mphi   ---- Mediator mass [eV]
        g      ---- Coupling. The interaction Lagrangian is
                    - (1/2) g  \bar{\psi}\psi \phi,
                    with \psi and \phi the neutrino and scalar fields,
                    respectively; and the 1/2 factor is for Majorana
                    neutrinos.
        mntot  ---- Sum of neutrino masses [eV]
        si     ---- Spectral index of the injected flux, given by
                    f(E) \propto E^{-si } exp(
        norm         ---- Normalization of the injected flux
                          This will be the free-streaming flavor-summed neutrino+antineutrino flux arriving at Earth at E=100 TeV
        """
        if mphi != None:
            self.my_object.mphi = mphi
        if g != None:
            self.my_object.g = g
        if mntot != None:
            self.my_object.mntot = mntot
        if si != None:
            self.my_object.si = si
        if norm != None:
            self.my_object.norm = norm
            
        self.evolved = False
        
    def evolve(self):
        """Evolve the neutrino flux"""
        self.evolved = True
        self.my_object.evolve()

    def get_flux(self):
        """Obtain the flux of each neutrino mass eigenstate, ordered from smallest to largest mass. Remember to call evolve() first"""
        N_bins = self.my_object.get_N_bins_E()
        flx = np.zeros([3, N_bins])
        if(not self.evolved):
            warnings.warn("You have not evolved the neutrino flux! Zero flux will be returned.")
            return flx

        for i in range(3):
            for j in range(N_bins):
                flx[i][j] = self.my_object.get_flux(i, j)
                
        return flx

    def get_flux_fla(self):
        """Obtain the flux of each neutrino flavour, ordered as {e, mu, tau}. Remember to call evolve() first"""
        N_bins = self.my_object.get_N_bins_E()
        flx = np.zeros([3, N_bins])
        if(not self.evolved):
            warnings.warn("You have not evolved the neutrino flux! Zero flux will be returned.")
            return flx
        
        for i in range(3):
            for j in range(N_bins):
                flx[i][j] = self.my_object.get_flux_fla(i, j)

        return flx

    def interp_flux_el(self, energy):
        """Obtain the nu_e flux at any given energy by interpolating"""
        return interp.interp1d(np.log10(self.get_energies()), self.get_flux_fla()[0] * self.get_energies()**self.my_object.si)(np.log10(energy)) / energy**self.my_object.si
    def interp_flux_mu(self, energy):
        """Obtain the nu_mu flux at any given energy by interpolating"""
        return interp.interp1d(np.log10(self.get_energies()), self.get_flux_fla()[1] * self.get_energies()**self.my_object.si)(np.log10(energy)) / energy**self.my_object.si
    def interp_flux_ta(self, energy):
        """Obtain the nu_tau flux at any given energy by interpolating"""
        return interp.interp1d(np.log10(self.get_energies()), self.get_flux_fla()[2] * self.get_energies()**self.my_object.si)(np.log10(energy)) / energy**self.my_object.si
            
    def get_energies(self):
        """Obtain the energy bin centers"""
        N_bins = self.my_object.get_N_bins_E()
        energies = np.zeros(N_bins)
        
        for i in range(N_bins):
            energies[i] = self.my_object.get_energy(i)

        return energies

    def check_energy_conservation(self):
        """Let E_FS and E_int be the total energy assuming free streaming
        and interactions, respectively. This function returns
           (E_int-E_FS)/E_FS"""
        return self.my_object.check_energy_conservation()
