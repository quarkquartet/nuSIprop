import nuSIprop

# Construct the object. Only the parameters mphi, g, mntot and si are mandatory
evolver = nuSIprop.pyprop(mphi = 6e5, # Mediator mass [eV]
			  g = 0.01, # Coupling
			  mntot = 0.1, # Sum of neutrino masses [eV]
			  si = 2.5, # Spectral index
			  norm = 5.991745031706099, # Normalization of the free-streaming flux at 100 TeV [Default = 1]
			  majorana = True, # Majorana neutrinos? [Default = True]
			  non_resonant = True, # Include non s-channel contributions? Relevant for couplings g>~0.1 [Default = True]
			  normal_ordering = True, # Normal neutrino mass ordering? [Default = True]
			  N_bins_E = 100, # Number of energy bins, uniformly distributed in log space [Default = 300]
			  lEmin = 9, # log_10 (E_min/eV) [Default = 13]
			  lEmax = 14, # log_10 (E_max/eV) [Default = 17]
			  zmax = 5, # Largest redshift at which sources are included [Default = 5]
			  flav = 2, # Flavor of interacting neutrinos [0=e, 1=mu, 2=tau. Default = 2]
			  phiphi = False # Consider double-scalar production? If set to true, the files xsec/alpha_phiphi.bin and xsec/alphatilde_phiphi.bin must exist [Default = False]
                          )

# Evolve it
evolver.evolve()
flx = evolver.get_flux_fla() # flx[0] is the nu_e flux, flx[1] the nu_mu flux, and flx[2] the nu_tau flux

# Output the result
print("#Energy[eV]  nu_e flux   nu_mu flux  nu_tau flux")
for (energy, flx_e, flx_mu, flx_ta) in zip(evolver.get_energies(),
                                           flx[0], flx[1], flx[2]):
    print("%.5e  %.4e  %.4e  %.4e" % (energy, flx_e, flx_mu, flx_ta))


# Other usage tips:

# -If you want the flux in the mass basis, just call
#  evolver.get_flux(). Then, flx[i] will be the flux of the neutrino
#  mass eigenstate with mass m[i]. Neutrino mass eigenstates are ordered
#  in the usual convention: nu_1 (nu_3) is the eigenstate with the smallest
#  (largest) electron neutrino content.

# -You can check that the total neutrino energy is conserved by
# checking the output of evolver.check_energy_conservation(). This
# returns the relative difference in total energy with and without
# neutrino self-interactions.  It thus gives an idea of the
# size of numerical errors

# -To change the physics parameters between consecutive runs, you
#  can call
#    evolver.set_parameters(mphi, g, mntot, si, norm)
#  to change respectively the mediator mass [eV], the coupling, the
#  total neutrino mass [eV], the spectral index of the propagated
#  flux, and its normalization.
#  Remember calling evolve() after updating them to obtain the correct
#  flux.

# -It is also possible to obtain the neutrino flux at an arbitrary
#  energy, by interpolating the evolved flux. For this, you can call
#    evolver.interp_flux_el(E)
#    evolver.interp_flux_mu(E)
#    evolver.interp_flux_ta(E)
#  to obtain the nu_e, nu_mu and nu_tau flux at the energy E [in eV]
