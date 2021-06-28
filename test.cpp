#include "nuSIprop.hpp"
#include <cstdio>

int main(){
  // Construct evolver object
  nuSIprop::calculate_flux evolver(6e5, // Mediator mass [eV]
				   0.01, // Coupling
				   0.1, // Sum of neutrino masses [eV]
				   2.5, // Spectral index
				   5.991745031706099, // Normalization of the free-streaming flux at 100 TeV [Default = 1]
				   true, // Majorana neutrinos? [Default = true]
				   true, // Include non s-channel contributions? Relevant for couplings g>~0.1 [Default = true]
				   true, // Normal neutrino mass ordering? [Default = true]
				   100, // Number of energy bins, uniformly distributed in log space [Default = 300]
				   9, // log_10 (E_min/eV) [Default = 13]
				   14, // log_10 (E_max/eV) [Default = 17]
				   5, // Largest redshift at which neutrino sources are included [Default = 5]
				   2, // Flavor of interacting neutrinos [0=e, 1=mu, 2=tau. Default = 2]
				   false // Consider double-scalar production? If set to true, the files xsec/alpha_phiphi.bin and xsec/alphatilde_phiphi.bin must exist [Default = false]
				   );

  // Evolve it
  evolver.evolve();

  // Output the result
  printf("#Energy[eV]  nu_e flux   nu_mu flux  nu_tau flux\n");
  for(int i=0; i<evolver.get_N_bins_E(); ++i)
    printf("%.5e  %.4e  %.4e  %.4e\n",
	   evolver.get_energy(i), // Energy bin centers
	   evolver.get_flux_fla(0, i), // nu_e flux
	   evolver.get_flux_fla(1, i), // nu_mu flux
	   evolver.get_flux_fla(2, i) // nu_tau flux
	   );
  
  /* Other usage tips:
     
     -If you want the flux in the mass basis, just call
      evolver.get_flux(i, j), where i={0,1,2} is the mass eigenstate
      and j is the energy bin. Neutrino mass eigenstates are ordered
      in the usual convention: nu_1 (nu_3) is the eigenstate with the smallest
      (largest) electron neutrino content.
     
     -You can check that the total neutrino energy is conserved by
     checking the output of evolver.check_energy_conservation(). This
     returns the relative difference in total energy with and without
     neutrino self-interactions.  It thus gives an idea of the
     size of numerical errors.

     -The public member variables
        evolver.mphi
	evolver.g
	evolver.mntot
	evolver.si
	evolver.norm
      can be changed between consecutive runs. They correspond to the
      mediator mass [eV], the coupling, the total neutrino mass [eV],
      the spectral index of the propagated flux, and its normalization.
      Remember calling evolve() after updating them to obtain the correct
      flux.
  */
}
