#ifndef _FLUX_H_
#define _FLUX_H_

#include "interp.hpp"
#include "aux.hpp"

#include <gsl/gsl_linalg.h>
#include <gsl/gsl_sf_dilog.h>
#include <complex.h>
#include <iostream>
#include <cmath>
#include <cstring>
#include "polylogarithm/src/cpp/Li2.hpp"
#include "polylogarithm/src/cpp/Li3.hpp"

#define SQR(x) ((x) * (x))       // square of a number
#define CUB(x) ((x) * (x) * (x)) // cube of a number

namespace nuSIprop
{

  class calculate_flux
  {
    /**
     * This class evolves an astrophysical neutrino flux assuming scalar
     * neutrino self-interactions with the cosmic neutrino background.
     *
     * The injected flux is assumed to follow a power law spectrum with a
     * high-energy cutoff, and a redshift dependence given by the Star
     * Formation Rate (see arXiv:0804.4008).
     *
     * Usage: after calling the evolve() method, get_flux(i, j) [0<i<3, 0<j<N_bins_E]
     *        will return the neutrino spectrum for each mass eigenstate at redshift 0
     *        and energy bin j.
     *        The function get_flux_fla(i, j) returns the neutrino spectrum
     *        in flavor space (i=0,1,2 corresponding to e,mu,tau)
     *
     * For convenience, we provide 3 additional functions:
     *        check_energy_conservation() : returns the relative difference in total energy
     *                                      with and without neutrino self-interactions.
     *                                      Since these conserve energy, this quantity should
     *                                      be small. It thus gives an idea of the numerical
     *                                      error
     *        get_energy(i): returns the [logarithmic] central energy of the bin i
     *        get_N_bins_E(): returns the number of energy bins
     *
     * As well as the public member variables:
     *        mphi  --- Mediator mass [eV]
     *        g     --- Yukawa coupling [eV]
     *        mntot --- Total neutrino mass [eV]
     *        si    --- Spectral index
     *        norm  --- Normalization of the free-streaming flux at 100 TeV
     * that can be modified between different runs
     *
     * Authors: Ivan Esteban, Sujata Pandey
     */

  public:
    calculate_flux() : calculate_flux(1e7, 0.1, 0.1, 2, 1, true, false) {} // We must declare a constructor without parameters for cython compatibility. For speed reasons, we set non_resonant to false. This will avoid loading the interpolating files

    calculate_flux(double mphi_, double g_, double mntot_, double si_,
                   double norm_ = 1,
                   bool majorana_ = true, bool non_resonant_ = true, bool normal_ordering_ = true,
                   int N_bins_E_ = 300, double lEmin_ = 12.0, double lEmax_ = 17.0,
                   double zmax_ = 5.0, int flav_ = 2, bool phiphi_ = false) : mphi(mphi_), g(g_), mntot(mntot_), si(si_),
                                                                              norm(norm_), majorana(majorana_), non_resonant(non_resonant_), normal_ordering(normal_ordering_),
                                                                              N_bins_E(N_bins_E_), lEmin(lEmin_), lEmax(lEmax_),
                                                                              flav(flav_), phiphi(phiphi_)
    {
      /**
       * Constructor
       *
       * Mandatory parameters:
       *   mphi_  ---- Mediator mass [eV]
       *   g_     ---- Coupling. The interaction Lagrangian is
       *                - (1/2) g  \bar{\psi}\psi \phi,
       *              with \psi and \phi the neutrino and scalar fields,
       *              respectively; and the 1/2 factor is for Majorana
       *              neutrinos.
       *   mntot_ ---- Sum of neutrino masses [eV]
       *   si_    ---- Spectral index of the injected flux, given by
       *
       * The energy dependence of the injected flux is given by
       *   f(E) = (E/100 TeV)^{-si_}
       *
       * Optional parameters:
       *   norm_           ---- Normalization of the injected flux. Default: 1
       *                        This will be the free-streaming flavor-summed neutrino+antineutrino flux arriving at Earth at E=100 TeV
       *   majorana        ---- Whether to consider Majorana (true) or Dirac (false) neutrinos. Default: true
       *   non_resonant    ---- Whether to consider non s- scattering channels (true) or not (false). Default: true
       *   normal_ordering ---- Whether to consider normal mass ordering (true) or inverted (false). Default: true
       *   N_bins_E_       ---- Number of energy bins. Default: 300
       *   lEmin_          ---- log_10(smallest energy considered / eV). Default: 12
       *   lEmax_          ---- log_10(largest energy considered / eV). Default: 17
       *   zmax_           ---- Largest redsfhit. Default: 5
       *   flav            ---- Flavor of interacting neutrinos: 0=e, 1=mu, 2=tau. Default: 2
       *   phiphi          ---- Whether to include double scalar production (true) or not (false).
       *                        This requires having the cross-section tables. Default: true
       *
       */
      // Allocate memory based on the number of energy bins
      flux = new double *[3]; // Neutrino spectrum in mass space
      for (int i = 0; i < 3; ++i)
        flux[i] = new double[N_bins_E];
      flux_fla = new double *[3]; // Neutrino spectrum in flavor space
      for (int i = 0; i < 3; ++i)
        flux_fla[i] = new double[N_bins_E];

      E_nu = new double[N_bins_E]; // Energy bin centers
      Emin = new double[N_bins_E]; // Smallest energy in each bin
      Emax = new double[N_bins_E]; // Largest energy in each bin

      for (int i = 0; i < N_bins_E; ++i)
      {
        Emin[i] = pow(10, lEmin + (lEmax - lEmin) * (i * 1.0) / N_bins_E);
        E_nu[i] = pow(10, lEmin + (lEmax - lEmin) * (i + 0.5) / N_bins_E);
        Emax[i] = pow(10, lEmin + (lEmax - lEmin) * (i + 1.0) / N_bins_E);
      }

      /* The redshift spacing is equal to the energy spacing.
       * I.e., 1+z[i] = (1+z[0]) * (Emax[j]/Emin[j])^i
       * This dramatically reduces the computation time
       */
      N_steps_z = log((1 + zmax_) / (1 + 0)) / log(Emax[0] / Emin[0]) + 2; // Number of redshift steps
      z = new double[N_steps_z];                                           // Redshift steps, ordered from z=0 to z=zmax
      for (int i = 0; i < N_steps_z; ++i)
        z[i] = (1 + 0) * pow(Emax[0] / Emin[0], i) - 1;
      zmax = z[N_steps_z - 1];

      /* Compute neutrino mixing*/
      double t12, t13, t23, dcp;
      if (normal_ordering)
      {
        t12 = 33.44 * (M_PI / 180); // theta12 in rad. Normal Ordering, NuFIT5.0
        t13 = 8.57 * (M_PI / 180);  // theta13 in rad. Normal Ordering, NuFIT5.0
        t23 = 49.0 * (M_PI / 180);  // theta23 in rad. Normal Ordering, NuFIT5.0
        dcp = 195.0 * (M_PI / 180); // delta_CP in rad. Normal Ordering, NuFIT5.0
      }
      else
      {
        t12 = 33.45 * (M_PI / 180); // theta12 in rad. Inverted Ordering, NuFIT5.0
        t13 = 8.61 * (M_PI / 180);  // theta13 in rad. Inverted Ordering, NuFIT5.0
        t23 = 49.3 * (M_PI / 180);  // theta23 in rad. Inverted Ordering, NuFIT5.0
        dcp = 286.0 * (M_PI / 180); // delta_CP in rad. Inverted Ordering, NuFIT5.0
      }
      // Sines and cosines of mixing angles
      double c12 = cos(t12);
      double c13 = cos(t13);
      double c23 = cos(t23);
      double s12 = sin(t12);
      double s13 = sin(t13);
      double s23 = sin(t23);
      std::complex<double> del(cos(dcp), sin(dcp));
      // Standard leptonic mixing matrix
      U[0][0] = c12 * c13;
      U[0][1] = s12 * c13;
      U[0][2] = s13 * 1.0 / del;
      U[1][0] = -s12 * c23 - c12 * s23 * s13 * del;
      U[1][1] = c12 * c23 - s12 * s23 * s13 * del;
      U[1][2] = s23 * c13;
      U[2][0] = s12 * s23 - c12 * c23 * s13 * del;
      U[2][1] = -c12 * s23 - s12 * c23 * s13 * del;
      U[2][2] = c23 * c13;

      // Set up the phi-phi cross section interpolators only if needed to save time
      if (non_resonant && phiphi)
      {
        spl_alphaTilde_phiphi = interp::spline_ND<2>({5000, 100}, "xsec/alphatilde_phiphi.bin", false, {true, false, false}, true);
        spl_alpha_phiphi = interp::spline_ND<3>({1000, 1000, 100}, "xsec/alpha_phiphi.bin", false, {true, false, false, false}, true);
      }
    }

    /* Parameters that can be modified after an object has been created (see documentation in the constructor for their meaning) */
    double mphi, g, mntot, si, norm;

    void evolve(void)
    {
      /**
       * Evolves the neutrino flux, updating the flux member array
       */

      /* Recompute quantities that may change from one run to another */
      // Neutrino masses
      double dmq21 = 7.42e-5; // [eV^2] Normal Ordering, NuFIT5.0
      double dmqAT;
      if (normal_ordering)
        dmqAT = 2.514e-3; // dmq31 [eV^2] Normal Ordering, NuFIT5.0
      else
        dmqAT = -2.497e-3; // dmq32 [eV^2] Normal Ordering, NuFIT5.0

      double mL = nuSIaux::getmL(mntot, dmq21, dmqAT);
      if (normal_ordering)
      {
        mn[0] = mL;
        mn[1] = sqrt(dmq21 + SQR(mL));
        mn[2] = sqrt(dmqAT + SQR(mL));
      }
      else
      {
        mn[2] = mL;
        mn[1] = sqrt(SQR(mL) - dmqAT);
        mn[0] = sqrt(SQR(mn[1]) - dmq21);
      }
      // Normalization
      norm_total = norm / flux_FS_E0();

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < N_bins_E; ++j)
          flux[i][j] = 0;

      // The evolution equations will be M*flux = v, and the solution will be stored in x
      gsl_matrix *M = gsl_matrix_alloc(3, 3);
      gsl_vector *v = gsl_vector_alloc(3);
      gsl_vector *x = gsl_vector_alloc(3);
      gsl_permutation *p = gsl_permutation_alloc(3); // Auxiliary variable needed by GSL

      // Auxiliary tables with integrals: we allocate them in the heap to avoid running out of memory
      double *tbl_Gamma = new double[N_bins_E + N_steps_z - 2];      // tbl_Gamma[i] = Gamma(Emin[i-k]*(1+z[k]), Emax[i-k]*(1+z[k]))
      double *tbl_alphaTilde = new double[N_bins_E + N_steps_z - 2]; // tbl_alphaTilde[i] = alphaTilde(Emin[i-k]*(1+z[k]), Emax[i-k]*(1+z[k]))
      double **tbl_alpha = new double *[N_bins_E + N_steps_z - 2];   // tbl_alpha[i][m] = alpha(Emin[i-k]*(1+z[k]), Emax[i-k]*(1+z[k]), Emin[m-k]*(1+z[k]), Emax[m-k]*(1+z[k]))
      for (int i = 0; i < N_bins_E + N_steps_z - 2; ++i)
      {
        double Emin_i, Emax_i;
        if (i < N_bins_E)
        {
          Emin_i = Emin[i];
          Emax_i = Emax[i];
        }
        else
        { // Auxiliary values for highest bins at highest redshifts
          Emin_i = Emin[N_bins_E - 1] * (1 + z[i - N_bins_E + 1]);
          Emax_i = Emax[N_bins_E - 1] * (1 + z[i - N_bins_E + 1]);
        }
        tbl_Gamma[i] = Gamma(Emin_i, Emax_i);
        tbl_alphaTilde[i] = alphaTilde(Emin_i, Emax_i);

        tbl_alpha[i] = new double[N_bins_E + N_steps_z - 2];
        for (int m = i + 1; m < N_bins_E + N_steps_z - 2; ++m)
        {
          double Emax_m, Emin_m;
          if (m < N_bins_E)
          {
            Emin_m = Emin[m];
            Emax_m = Emax[m];
          }
          else
          { // Auxiliary values for highest bins at highest redshifts
            Emin_m = Emin[N_bins_E - 1] * (1 + z[m - N_bins_E + 1]);
            Emax_m = Emax[N_bins_E - 1] * (1 + z[m - N_bins_E + 1]);
          }
          tbl_alpha[i][m] = alpha(Emin_i, Emax_i, Emin_m, Emax_m);
        }
      }

      double alpha_wo_mixing[N_bins_E];
      double dlogz = log(1 + z[1]) - log(1 + z[0]);
      for (int i = N_steps_z - 1; i > 0; --i) // Loop for redshift, starting from z=zmax to z=0. We will obtain the solution of the evolution equations at z[i-1]
      {
        double H = get_H(z[i - 1]); // Hubble parameter

        double alpha_cum[3] = {0, 0, 0};
        // If we are only considering resonant contributions, alpha \propto tminus-tplus.
        // Thus, when we sum over alpha we are computing many times essentially the same quantity.
        // To speed everything up, the variable alpha_cum keeps track of alpha

        for (int j = N_bins_E; j > 0; --j)
        {                                                                                                 // Loop over energies, starting from highest energies. Current energy bin = E_nu[j-1]
          double Gamma_wo_mixing = get_nd(z[i - 1]) / SQR(1 + z[i - 1]) * tbl_Gamma[j + i - 2];           // Gamma(Emin[j-1]*(1+z[i-1]), Emax[j-1]*(1+z[i-1]));
          double alphaTilde_wo_mixing = get_nd(z[i - 1]) / SQR(1 + z[i - 1]) * tbl_alphaTilde[j + i - 2]; // alphaTilde(Emin[j-1]*(1+z[i-1]), Emax[j-1]*(1+z[i-1]));
          if (non_resonant)
            for (int m = j; m < N_bins_E; ++m)
              alpha_wo_mixing[m] = get_nd(z[i - 1]) / SQR(1 + z[i - 1]) * tbl_alpha[j + i - 2][m + i - 1]; // alpha(Emin[j-1]*(1+z[i-1]), Emax[j-1]*(1+z[i-1]), Emin[m]*(1+z[i-1]), Emax[m]*(1+z[i-1]));
          else if (j != N_bins_E)
          {
            alpha_wo_mixing[j] = get_nd(z[i - 1]) / SQR(1 + z[i - 1]) * tbl_alpha[j + i - 2][j + i - 1]; // alpha(Emin[j-1]*(1+z[i-1]), Emax[j-1]*(1+z[i-1]), Emin[j]*(1+z[i-1]), Emax[j]*(1+z[i-1]));
            for (int l = 0; l < 3; ++l)                                                                  // Loop over flavors
              alpha_cum[l] += flux[l][j] * alpha_wo_mixing[j] / (Emax[j] - Emin[j]) / (Emax[j - 1] - Emin[j - 1]);
          }

          /* We fill the matrices representing the evolution equations */
          for (int k = 0; k < 3; ++k)
          {
            double src = (1 + z[i - 1]) * dlogz / H * Lum(z[i], Emin[j - 1], Emax[j - 1], k); // Source term

            if (!non_resonant && j != N_bins_E)
              for (int l = 0; l < 3; ++l) // Loop over flavors
                src += (1 + z[i - 1]) * dlogz / H * alpha_cum[l] * std::norm(U[flav][k]) * std::norm(U[flav][l]) * (Emax[j - 1] - Emin[j - 1]);
            else
              for (int m = j; m < N_bins_E; ++m) // Loop over bins with energies er[m], with m>=j
                for (int l = 0; l < 3; ++l)      // Loop over flavors
                  src += (1 + z[i - 1]) * dlogz / H * flux[l][m] * alpha_wo_mixing[m] * std::norm(U[flav][k]) * std::norm(U[flav][l]) / (Emax[m] - Emin[m]);

            double Znr = flux[k][j - 1] + src;                                                                                                                                           // Numerator of the rhs of the evolution equation
            double Zdr = 1.0 + (1 + z[i - 1]) * dlogz / H * (Gamma_wo_mixing * std::norm(U[flav][k]) - alphaTilde_wo_mixing * SQR(std::norm(U[flav][k]))) / (Emax[j - 1] - Emin[j - 1]); // Denominator of the rhs of the evolution equation

            gsl_vector_set(v, k, Znr / Zdr);
            for (int l = 0; l < 3; ++l)
            {
              if (k == l)
                gsl_matrix_set(M, k, l, 1.0);
              else
                gsl_matrix_set(M, k, l,
                               (alphaTilde_wo_mixing * std::norm(U[flav][k]) * std::norm(U[flav][l]) / (Emax[j - 1] - Emin[j - 1])) / Zdr);
            }
          }

          /* Solve the evolution equations */
          int s;
          gsl_linalg_LU_decomp(M, p, &s);
          gsl_linalg_LU_solve(M, p, v, x);

          for (int k = 0; k < 3; ++k)
            flux[k][j - 1] = gsl_vector_get(x, k);
        }
      }

      // Free the allocated memory
      gsl_matrix_free(M);
      gsl_vector_free(v);
      gsl_vector_free(x);
      gsl_permutation_free(p);
      delete[] tbl_Gamma;
      delete[] tbl_alphaTilde;
      for (int i = 0; i < N_bins_E + N_steps_z - 2; ++i)
        delete[] tbl_alpha[i];
      delete[] tbl_alpha;

      // Divide by energy bin size
      for (int i = 0; i < N_bins_E; ++i)
        for (int k = 0; k < 3; ++k)
          flux[k][i] /= (Emax[i] - Emin[i]);

      // Convert to flavor space
      for (int i = 0; i < N_bins_E; ++i)
        for (int k = 0; k < 3; ++k)
          flux_fla[k][i] = std::norm(U[k][0]) * flux[0][i] + std::norm(U[k][1]) * flux[1][i] + std::norm(U[k][2]) * flux[2][i];
    }

    double check_energy_conservation(void)
    {
      /**
       * Let E_FS and E_int be the total energy assuming free streaming
       * and interactions, respectively. This function returns
       *    (E_int-E_FS)/E_FS
       */

      double E_FS = energy_FS();

      evolve();
      double E_int = 0;
      // \int dE * E * flux(E) = \int d(logE) * E^2 * flux(E)
      for (int i = 0; i < N_bins_E; ++i)
        for (int k = 0; k < 3; ++k)
          E_int += (log(Emax[i]) - log(Emin[i])) * SQR(E_nu[i]) * flux[k][i];

      return (E_int - E_FS) / E_FS;
    }

    inline double get_flux(int i, int j)
    {
      /**
       * Returns the neutrino flux for the mass eigenstate i and energy j
       */
      if ((i < 0) || (i >= 3))
      {
        std::cerr << "You asked for the flux of the mass eigenstate " << i << ", not in [0,1,2]. Zero will be returned." << std::endl;
        return 0;
      }
      else if (j < 0)
      {
        std::cerr << "You asked for the flux at the energy bin " << j << "<0! Zero will be returned." << std::endl;
        return 0;
      }
      else if (j > N_bins_E)
      {
        std::cerr << "You asked for the flux at the energy bin " << j << ", but there are only " << N_bins_E << " bins! Zero will be returned." << std::endl;
        return 0;
      }

      return flux[i][j];
    }

    inline double get_flux_fla(int i, int j)
    {
      /**
       * Returns the neutrino flux for the flavor i and energy j
       */
      if ((i < 0) || (i >= 3))
      {
        std::cerr << "You asked for the flux of the flavor eigenstate " << i << ", not in [0,1,2]. Zero will be returned." << std::endl;
        return 0;
      }
      else if (j < 0)
      {
        std::cerr << "You asked for the flux at the energy bin " << j << "<0! Zero will be returned." << std::endl;
        return 0;
      }
      else if (j > N_bins_E)
      {
        std::cerr << "You asked for the flux at the energy bin " << j << ", but there are only " << N_bins_E << " bins! Zero will be returned." << std::endl;
        return 0;
      }

      return flux_fla[i][j];
    }

    inline int get_N_bins_E()
    {
      return N_bins_E;
    }

    inline double get_energy(int i)
    {
      /**
       * Returns the [logarithmic] central energy of the bin i
       */
      if (i < 0)
      {
        std::cerr << "You asked for the energy at the bin " << i << "<0! Zero will be returned." << std::endl;
        return 0;
      }
      else if (i > N_bins_E)
      {
        std::cerr << "You asked for the energy at the bin " << i << ", but there are only " << N_bins_E << " bins! Zero will be returned." << std::endl;
        return 0;
      }

      return E_nu[i];
    }

    /* Constructors and destructors that take proper care of the dynamic memory */

    // Copy constructor
    calculate_flux(const calculate_flux &orig) : mphi(orig.mphi), g(orig.g), mntot(orig.mntot), si(orig.si),
                                                 norm(orig.norm),
                                                 majorana(orig.majorana), non_resonant(orig.non_resonant), normal_ordering(orig.normal_ordering),
                                                 N_bins_E(orig.N_bins_E), N_steps_z(orig.N_steps_z), lEmin(orig.lEmin), lEmax(orig.lEmax),
                                                 zmax(orig.zmax), flav(orig.flav), phiphi(orig.phiphi),
                                                 spl_alpha_phiphi(orig.spl_alpha_phiphi), spl_alphaTilde_phiphi(orig.spl_alphaTilde_phiphi)
    {
      // Properly take care of the memory in the arrays
      flux = new double *[3];
      for (int i = 0; i < 3; ++i)
        flux[i] = new double[N_bins_E];
      flux_fla = new double *[3];
      for (int i = 0; i < 3; ++i)
        flux_fla[i] = new double[N_bins_E];
      E_nu = new double[N_bins_E];
      Emin = new double[N_bins_E];
      Emax = new double[N_bins_E];
      z = new double[N_steps_z];

      for (int i = 0; i < 3; ++i)
        std::memcpy(flux[i], orig.flux[i], N_bins_E * sizeof(double));
      for (int i = 0; i < 3; ++i)
        std::memcpy(flux_fla[i], orig.flux_fla[i], N_bins_E * sizeof(double));
      std::memcpy(E_nu, orig.E_nu, N_bins_E * sizeof(*orig.E_nu));
      std::memcpy(Emin, orig.Emin, N_bins_E * sizeof(*orig.Emin));
      std::memcpy(Emax, orig.Emax, N_bins_E * sizeof(*orig.Emax));
      std::memcpy(z, orig.z, N_steps_z * sizeof(*orig.z));

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          U[i][j] = orig.U[i][j];
    }

    // Assignment operator
    calculate_flux &operator=(const calculate_flux &rhs)
    {
      mphi = rhs.mphi;
      g = rhs.g;
      mntot = rhs.mntot;
      si = rhs.si;
      norm = rhs.norm;
      majorana = rhs.majorana;
      non_resonant = rhs.non_resonant;
      normal_ordering = rhs.normal_ordering;
      N_bins_E = rhs.N_bins_E;
      N_steps_z = rhs.N_steps_z;
      lEmin = rhs.lEmin;
      lEmax = rhs.lEmax;
      zmax = rhs.zmax;
      flav = rhs.flav;
      phiphi = rhs.phiphi;
      spl_alpha_phiphi = rhs.spl_alpha_phiphi;
      spl_alphaTilde_phiphi = rhs.spl_alphaTilde_phiphi;

      // Properly take care of the memory in the arrays
      for (int i = 0; i < 3; ++i)
        delete[] flux[i];
      delete[] flux;
      for (int i = 0; i < 3; ++i)
        delete[] flux_fla[i];
      delete[] flux_fla;
      delete[] E_nu;
      delete[] Emin;
      delete[] Emax;
      delete[] z;

      flux = new double *[3];
      for (int i = 0; i < 3; ++i)
        flux[i] = new double[N_bins_E];
      flux_fla = new double *[3];
      for (int i = 0; i < 3; ++i)
        flux_fla[i] = new double[N_bins_E];
      E_nu = new double[N_bins_E];
      Emin = new double[N_bins_E];
      Emax = new double[N_bins_E];
      z = new double[N_steps_z];

      for (int i = 0; i < 3; ++i)
        std::memcpy(flux[i], rhs.flux[i], N_bins_E * sizeof(double));
      for (int i = 0; i < 3; ++i)
        std::memcpy(flux_fla[i], rhs.flux_fla[i], N_bins_E * sizeof(double));
      std::memcpy(E_nu, rhs.E_nu, N_bins_E * sizeof(*rhs.E_nu));
      std::memcpy(Emin, rhs.Emin, N_bins_E * sizeof(*rhs.Emin));
      std::memcpy(Emax, rhs.Emax, N_bins_E * sizeof(*rhs.Emax));
      std::memcpy(z, rhs.z, N_steps_z * sizeof(*rhs.z));

      for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
          U[i][j] = rhs.U[i][j];

      return *this;
    }

    // Destructor
    ~calculate_flux()
    {
      for (int i = 0; i < 3; ++i)
        delete[] flux[i];
      delete[] flux;
      for (int i = 0; i < 3; ++i)
        delete[] flux_fla[i];
      delete[] flux_fla;
      delete[] E_nu;
      delete[] Emin;
      delete[] Emax;
      delete[] z;
    }

  private:
    /* Constant quantities, not to be modified by the external user between two runs*/
    double **flux;                // Matrix that will contain the final neutrino spectrum in mass space
    double **flux_fla;            // Matrix that will contain the final neutrino spectrum in flavor space
    double *E_nu;                 // Array that will contain the energy bin centers
    double *Emin, *Emax;          // Smallest and largest energies in each bin
    double *z;                    // Redshift bins
    double E0 = 1e14;             // Pivot energy in the flux [eV]
    int N_integ_z = 100;          // Number of redshift bins to integrate the free-streaming energy
    std::complex<double> U[3][3]; // Leptonic mixing matrix
    // See the documentation of the parameters below in the constructor
    bool majorana, non_resonant, normal_ordering;
    int N_bins_E;
    int N_steps_z;
    double lEmin, lEmax;
    double zmax;
    int flav;
    bool phiphi;

    // Interpolators
    interp::spline_ND<3> spl_alpha_phiphi;
    interp::spline_ND<2> spl_alphaTilde_phiphi;

    /* Parameters that are derived from the user-modifiable parameters
       These are modified each time evolve() is called, so there's no need to take care of them in constructors
       */
    double mn[3];      // Individual neutrino masses
    double norm_total; // Parameter encoding the global normalization of the flux

    /* Physics functions */

    static inline double get_nd(double z)
    {
      /**
       * Returns the cosmic *neutrino* number density **of each mass eigenstate** at redshift z [eV^3]
       * To get the total flavor-summed neutrino+antineutrino number density, multiply this by 6
       */
      return 4.3528e-13 * CUB(1 + z);
    }

    static inline double get_H(double z)
    {
      /**
       * Returns the Hubble parameter at redshift z [eV]
       * We take H0 = 70 km/s/Mpc, Omega_M = 0.308, Omega_L = 0.692
       */
      return 1.5e-33 * pow(0.692 + 0.308 * CUB(1 + z), 0.5);
    }

    static inline double get_SFR(double z)
    {
      /**
       * Returns the Star Formation Rate at redshift z,
       * using the parametrization in
       * Yüksel, Kistler, Beacom & Hopkins, arXiv:0804.4008
       *
       * The normalization is arbitrary
       */

      return pow(pow(1 + z, -3.4 * 10) +
                     pow((1 + z) / 5161, 0.3 * 10) +
                     pow((1 + z) / 9.06, 3.5 * 10),
                 -1. / 10.);
    }

    static inline double RSN(double z)
    {
      /**
       * Returns the R_SN(z)
       * See eq 4 of 0812.3157, or eq 6 of 1004.3311
       * The R_SF(z) or \dot{\rho}_* (z) refers to the SFR defined above.
       */
      double m_solar = 1.989 * 56.1; // Solar mass in unit 1e64 eV.
      return get_SFR(z) / (0.01 * m_solar);
    }

    static inline double dNdE(double E)
    {
      /**
       * Returns the energy spectrum of Fermi-Dirac distribution.
       */
      double Etot = 3 * 6.24; // E_nu^tot in the unit 1e64 eV.
      double Tnue = 6e6;      // 6 MeV of DSNB temperature
      return Etot * 120 * pow(E, 2) / (6 * 7 * pow(M_PI, 4) * pow(Tnue, 4) * (exp(E / Tnue) + 1));
    }

    static inline double Li_2(double x)
    {
      return polylogarithm::Li2(x);
    }

    static inline double Li_3(double x)
    {
      return polylogarithm::Li3(x);
    }

    static inline double Lum_int(double z, double E)
    {
      /**
       * Return the integration of Fermi-Dirac distribution.
       *  */
      double Etot = 3 * 6.24; // E_nu^tot in the unit 1e64 eV.
      double Tnue = 6e6;
      return (Etot * 120 / (6 * 7 * pow(M_PI, 4) * pow(Tnue, 2))) * (-E * E * (1 + z) * log(exp(-E * (1 + z) / Tnue) + 1) / Tnue + 2 * E * Li_2(-exp(-E * (1 + z) / Tnue)) + 2 * Tnue * Li_3(-exp(-E * (1 + z) / Tnue)) / (1 + z));
    }

    // inline double Lum(double z, double Em, double Ep, int i)
    // {
    //   /**
    //    * Returns
    //    *   \int_Em^Ep L^i(z, E*(1+z)) dE
    //    * where L^i(z, E*(1+z)) is the source distribution function for the mass eigenstate i=0,1,2
    //    */

    //   return norm_total / 3.0 * get_SFR(z) * (Ep * pow(Ep / E0 * (1 + z), -si) - Em * pow(Em / E0 * (1 + z), -si)) / (1 - si);
    // }

    inline double Lum(double z, double Em, double Ep, int i)
    {
      return (Lum_int(z, Ep) - Lum_int(z, Em)) * RSN(z);
    }

    /* Auxiliary functions to compute fluxes and energies */

    double flux_FS_E0()
    {
      /**
       * Returns the free streaming flux at Earth at the pivot energy E0
       * This is given by the integral of
       *    \int_0^zmax dz/H(z) \sum_i L^i(z, E0*(1+z))
       * where L^i(z, E*(1+z)) is the source distribution function.
       */

      double res = 0;
      double z_min = 0, z_max = zmax;

      for (int f = 0; f < N_integ_z; f++)
      {
        // Limits of the integral
        double a = z_min + f * (z_max - z_min) / N_integ_z;
        double b = z_min + (f + 1.0) * (z_max - z_min) / N_integ_z;

        // Nodes at which the integrand will be evaluated
        double z1 = (b - a) / 2. * nuSIaux::x_integ[0] + (b + a) / 2.,
               z2 = (b - a) / 2. * nuSIaux::x_integ[1] + (b + a) / 2.,
               z3 = (b - a) / 2. * nuSIaux::x_integ[2] + (b + a) / 2.;

        res += (b - a) / 2. * (nuSIaux::w_integ[0] * pow(1 + z1, -si) * get_SFR(z1) / get_H(z1) + nuSIaux::w_integ[1] * pow(1 + z2, -si) * get_SFR(z2) / get_H(z2) + nuSIaux::w_integ[2] * pow(1 + z3, -si) * get_SFR(z3) / get_H(z3));
      }
      return res;
    }

    double energy_FS()
    {
      /**
       * Returns the total energy for free streaming flux, propagating neutrinos from z=zmax to z=0
       * This is given by the integral of
       *    \int_0^zmax dz/H(z) \int_Emin^Emax dE E \sum_i L^i(z, E*(1+z))
       * where L^i(z, E*(1+z)) is the source distribution function.
       */
      return energy_FS(0, zmax);
    }

    double energy_FS(double z_min, double z_max)
    {
      /**
       * Returns the total energy for free streaming flux, propagating neutrinos from z_max to z_min
       * This is given by the integral of
       *    \int_zmin^zmax dz/H(z) \int_Emin^Emax dE E \sum_i L^i(z, E*(1+z))
       * where L^i(z, E*(1+z)) is the source distribution function.
       */

      double res = 0;
      for (int f = 0; f < N_integ_z; f++)
      {
        // Limits of the integral
        double a = z_min + f * (z_max - z_min) / N_integ_z;
        double b = z_min + (f + 1.0) * (z_max - z_min) / N_integ_z;

        // Nodes at which the integrand will be evaluated
        double z1 = (b - a) / 2. * nuSIaux::x_integ[0] + (b + a) / 2.,
               z2 = (b - a) / 2. * nuSIaux::x_integ[1] + (b + a) / 2.,
               z3 = (b - a) / 2. * nuSIaux::x_integ[2] + (b + a) / 2.;

        res += (b - a) / 2. * (nuSIaux::w_integ[0] * Lum_times_E(z1, pow(10, lEmin), pow(10, lEmax)) / get_H(z1) + nuSIaux::w_integ[1] * Lum_times_E(z2, pow(10, lEmin), pow(10, lEmax)) / get_H(z2) + nuSIaux::w_integ[2] * Lum_times_E(z3, pow(10, lEmin), pow(10, lEmax)) / get_H(z3));
      }
      return res;
    }

    double Lum_times_E(double z, double Em, double Ep)
    {
      /**
       * Returns
       *   \int_Em^Ep E \sum_i L(z, E*(1+z), i) dE
       * where L(z, E*(1+z), i) is the source distribution function.
       * We will assume that all flavors are emitted equally
       */

      if (fabs(si - 2) < 1e-5) // Possible roundoff uncertainties: Taylor-expand
        return norm_total * get_SFR(z) * pow(E0 / (1 + z), si) * (log(Ep / Em) + (2 - si) / 2.0 * (SQR(log(Ep)) - SQR(log(Em))));
      else
        return norm_total * get_SFR(z) * pow(E0 / (1 + z), si) * (pow(Ep, 2 - si) - pow(Em, 2 - si)) / (2 - si);
    }

    /* Quantities related to the interaction */

    double scalar_width()
    {
      /**
       * Returns the scalar decay width
       */
      if (majorana)
        return SQR(g) * mphi / (16.0 * M_PI);
      else
        return SQR(g) * mphi / (8.0 * M_PI);
    }

    double Gamma(double Em, double Ep)
    {
      /**
       * Returns
       *  \sum_j \int_Em^Ep sigma_{ij}(E)/|U_{flav i}|^2 dE
       * where
       *  sigma_{ij} is the absorption cross section for an initial neutrino with mass i
       *   hitting on a cosmic neutrino with mass j
       *  i, j = 0,1,2
       */

      double Ga = scalar_width(); // Scalar decay width

      double tot = 0;
      for (int j = 0; j < 3; ++j)
      {
        /* Integration limits in terms of s/mphi^2 */
        double splus = 2 * mn[j] * Ep / SQR(mphi);
        double sminus = 2 * mn[j] * Em / SQR(mphi);

        /* s-channel */
        double Gamma_s;
        if (splus < 1e-5) // We Taylor-expand atandiff to avoid roundoff errors
          Gamma_s = SQR(SQR(g)) / (32 * M_PI * SQR(mphi) * Ga) *
                    (2 * mphi * ((Ga / mphi * (1 + SQR(Ga / mphi) + 2 * sminus)) / SQR(1 + SQR(Ga / mphi)) * (splus - sminus) + (Ga / mphi) / SQR(1 + SQR(Ga / mphi)) * SQR(splus - sminus)) + Ga * (log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * splus * (splus - 2)) -
                                                                                                                                                                                                     log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * sminus * (sminus - 2))));
        else
          Gamma_s = SQR(SQR(g)) / (32 * M_PI * SQR(mphi) * Ga) *
                    (2 * mphi * nuSIaux::atandiff(mphi * (splus - 1) / Ga, mphi * (sminus - 1) / Ga) + Ga * (log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * splus * (splus - 2)) -
                                                                                                             log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * sminus * (sminus - 2))));
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 \sum_{kl} |U_{flav k}|^2 |U_{flav l}|^2
        Gamma_s *= std::norm(U[flav][j]);
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_s;

        if (!non_resonant)
          continue;

        /* t-channel + u-channel (without interference) */
        double Gamma_t_u = SQR(SQR(g)) / (16 * M_PI * SQR(mphi)) *
                           (2 * log1p(splus) / splus - 2 * log1p(sminus) / sminus + log1p(splus) - log1p(sminus));
        if (Gamma_t_u < 0)
        { // Roundoff errors! Compute the integral numerically
          // Limits of the integral
          double a = sminus;
          double b = splus;
          // Nodes at which the integrand will be evaluated
          double z1 = (b - a) / 2. * nuSIaux::x_integ[0] + (b + a) / 2.,
                 z2 = (b - a) / 2. * nuSIaux::x_integ[1] + (b + a) / 2.,
                 z3 = (b - a) / 2. * nuSIaux::x_integ[2] + (b + a) / 2.;
          Gamma_t_u = SQR(SQR(g)) / (16 * M_PI * SQR(mphi)) *
                      (b - a) / 2. * (nuSIaux::w_integ[0] * ((z1 + 2) / (z1 * (z1 + 1)) - 2 / SQR(z1) * log1p(z1)) + nuSIaux::w_integ[1] * ((z2 + 2) / (z2 * (z2 + 1)) - 2 / SQR(z2) * log1p(z2)) + nuSIaux::w_integ[2] * ((z3 + 2) / (z3 * (z3 + 1)) - 2 / SQR(z3) * log1p(z3)));
        }
        // Prefactor: 2 * |U_{flav i}|^2 |U_{flav j}|^2 \sum_{kl} |U_{flav k}|^2 |U_{flav l}|^2
        if (majorana)
          Gamma_t_u *= 2 * std::norm(U[flav][j]);
        else
          Gamma_t_u *= 2 * std::norm(U[flav][j]);
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_t_u;

        /* t-u interference */
        double Gamma_tu = SQR(SQR(g)) / (32 * M_PI * SQR(mphi) * sminus * splus) *
                          (sminus * log1p(splus) * (2 + 2 * splus + splus * log(2 + splus)) -
                           splus * log1p(sminus) * (2 + 2 * sminus + sminus * log(2 + sminus)) +
                           sminus * splus * (nuSIaux::dilog1mdiff(splus, sminus) + nuSIaux::dilogdiff(splus, sminus)));
        if (Gamma_tu < 0)
        { // Roundoff errors! Compute the integral numerically
          // Limits of the integral
          double a = sminus;
          double b = splus;
          // Nodes at which the integrand will be evaluated
          double z1 = (b - a) / 2. * nuSIaux::x_integ[0] + (b + a) / 2.,
                 z2 = (b - a) / 2. * nuSIaux::x_integ[1] + (b + a) / 2.,
                 z3 = (b - a) / 2. * nuSIaux::x_integ[2] + (b + a) / 2.;
          Gamma_tu = SQR(SQR(g)) / (16 * M_PI * SQR(mphi)) *
                     (b - a) / 2. * (nuSIaux::w_integ[0] * (1 / z1 - 2 * (1 + z1) / (SQR(z1) * (2 + z1)) * log1p(z1)) + nuSIaux::w_integ[1] * (1 / z2 - 2 * (1 + z2) / (SQR(z2) * (2 + z2)) * log1p(z2)) + nuSIaux::w_integ[2] * (1 / z3 - 2 * (1 + z3) / (SQR(z3) * (2 + z3)) * log1p(z3)));
        }
        // Prefactor: (1/2) * |U_{flav i}|^2 |U_{flav j}|^2 \sum_{kl} |U_{flav k}|^2 |U_{flav l}|^2 . The 1/2 factor is for Dirac, that has half the amount of targets in the u-channel
        if (majorana)
          Gamma_tu *= std::norm(U[flav][j]);
        else
          Gamma_tu *= 0.5 * std::norm(U[flav][j]);
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_tu;

        /* s-t interference */
        double ga_red = Ga / mphi; // Reduced decay width
        /* Compute dilogarithms */
        // Auxiliary complex numbers
        double _Complex z1_plus = I * (1 + splus) / (2 * I + ga_red);
        double _Complex z1_minus = I * (1 + sminus) / (2 * I + ga_red);
        double _Complex z2_plus = conj(z1_plus);
        double _Complex z2_minus = conj(z1_minus);
        double _Complex dilogdiff_z1;
        double _Complex dilogdiff_z2;

        if (splus < 1e-5)
        { // We Taylor-expand dilogdiff to be fast and avoid roundoff errors
          dilogdiff_z1 = SQR(sminus) * (-I / 2 / (I + ga_red) - clog((I + ga_red) / (2 * I + ga_red)) / 2.) +
                         sminus * clog((I + ga_red) / (2 * I + ga_red)) - splus * clog((I + ga_red) / (2 * I + ga_red)) +
                         (SQR(splus) * (I / (I + ga_red) + clog((I + ga_red) / (2 * I + ga_red)))) / 2.;
          dilogdiff_z2 = SQR(sminus) * (I / 2 / (-I + ga_red) - clog((-I + ga_red) / (-2 * I + ga_red)) / 2.) +
                         sminus * clog((-I + ga_red) / (-2 * I + ga_red)) - splus * clog((-I + ga_red) / (-2 * I + ga_red)) +
                         (SQR(splus) * (-I / (-I + ga_red) + clog((-I + ga_red) / (-2 * I + ga_red)))) / 2.;
        }
        else
        {
          dilogdiff_z1 = nuSIaux::dilogdiff_complex(z1_plus, z1_minus);
          dilogdiff_z2 = nuSIaux::dilogdiff_complex(z2_plus, z2_minus);
        }

        double Gamma_st = -SQR(SQR(g)) / (32 * M_PI * SQR(mphi) * (1 + SQR(ga_red))) *
                          (creal(dilogdiff_z1) + creal(dilogdiff_z2) + ga_red * (cimag(dilogdiff_z2) - cimag(dilogdiff_z1)) + 2 * ga_red * carg(1 - z2_plus) * log1p(splus) - 2 * ga_red * carg(1 - z2_minus) * log1p(sminus) + log1p(4 / SQR(ga_red)) * (log1p(sminus) - log1p(splus)) + log1p(SQR(-1 + splus) / SQR(ga_red)) * log1p(splus) - log1p(SQR(-1 + sminus) / SQR(ga_red)) * log1p(sminus) + (1 + SQR(ga_red)) * (log1p(SQR(-1 + sminus) / SQR(ga_red)) - log1p(SQR(-1 + splus) / SQR(ga_red))) + 2 * nuSIaux::dilogdiff(splus, sminus));
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 \sum_{kl} |U_{flav k}|^2 |U_{flav l}|^2
        Gamma_st *= std::norm(U[flav][j]);
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_st;

        /* s-u interference (only for Majorana fermions) */
        double Gamma_su = 0;
        if (majorana)
          Gamma_su = Gamma_st;
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_su;

        /* Double scalar production */
        double Gamma_pp = 0;
        if ((splus > 4) && phiphi)
        {
          if (sminus > 4)
            Gamma_pp = SQR(SQR(g)) / (128. * M_PI * SQR(mphi)) * (12 * sqrt((-4 + sminus) / sminus) - 12 * sqrt((-4 + splus) / splus) - 2 * log(SQR(sqrt(-4 + sminus) - sqrt(sminus)) / 4.) * log(SQR(-2 + sminus + sqrt((-4 + sminus) * sminus)) / 4.) - ((6 + sminus * log((-2 + sminus) * sminus)) * log(SQR(-2 + sminus + sqrt((-4 + sminus) * sminus)) / SQR(2 - sminus + sqrt((-4 + sminus) * sminus)))) / sminus - 24 * (sqrt((-4 + sminus) / sminus) - sqrt((-4 + splus) / splus) - log(sqrt(-4 + sminus) + sqrt(sminus)) + log(sqrt(-4 + splus) + sqrt(splus))) + 2 * log(SQR(sqrt(-4 + splus) - sqrt(splus)) / 4.) * log(SQR(-2 + splus + sqrt((-4 + splus) * splus)) / 4.) + ((6 + splus * log((-2 + splus) * splus)) * log(SQR(-2 + splus + sqrt((-4 + splus) * splus)) / SQR(2 - splus + sqrt((-4 + splus) * splus)))) / splus + 8 * nuSIaux::dilogdiff(4 / SQR(sqrt(-4 + sminus) + sqrt(sminus)), 4 / SQR(sqrt(-4 + splus) + sqrt(splus))) + 2 * nuSIaux::dilogdiff(4 / SQR(-2 + sminus + sqrt((-4 + sminus) * sminus)), 4 / SQR(-2 + splus + sqrt((-4 + splus) * splus))));
          else
            Gamma_pp = SQR(SQR(g)) / (128. * M_PI * SQR(mphi)) * (12 * sqrt((-4 + 4) / 4) - 12 * sqrt((-4 + splus) / splus) - 2 * log(SQR(sqrt(-4 + 4) - sqrt(4)) / 4.) * log(SQR(-2 + 4 + sqrt((-4 + 4) * 4)) / 4.) - ((6 + 4 * log((-2 + 4) * 4)) * log(SQR(-2 + 4 + sqrt((-4 + 4) * 4)) / SQR(2 - 4 + sqrt((-4 + 4) * 4)))) / 4 - 24 * (sqrt((-4 + 4) / 4) - sqrt((-4 + splus) / splus) - log(sqrt(-4 + 4) + sqrt(4)) + log(sqrt(-4 + splus) + sqrt(splus))) + 2 * log(SQR(sqrt(-4 + splus) - sqrt(splus)) / 4.) * log(SQR(-2 + splus + sqrt((-4 + splus) * splus)) / 4.) + ((6 + splus * log((-2 + splus) * splus)) * log(SQR(-2 + splus + sqrt((-4 + splus) * splus)) / SQR(2 - splus + sqrt((-4 + splus) * splus)))) / splus + 8 * nuSIaux::dilogdiff(4 / SQR(sqrt(-4 + 4) + sqrt(4)), 4 / SQR(sqrt(-4 + splus) + sqrt(splus))) + 2 * nuSIaux::dilogdiff(4 / SQR(-2 + 4 + sqrt((-4 + 4) * 4)), 4 / SQR(-2 + splus + sqrt((-4 + splus) * splus))));

          if (Gamma_pp < 0)
          { // Roundoff errors! Compute the integral numerically
            // Limits of the integral
            double a = (sminus < 4) ? 4 : sminus;
            double b = splus;
            // Nodes at which the integrand will be evaluated
            double z1 = (b - a) / 2. * nuSIaux::x_integ[0] + (b + a) / 2.,
                   z2 = (b - a) / 2. * nuSIaux::x_integ[1] + (b + a) / 2.,
                   z3 = (b - a) / 2. * nuSIaux::x_integ[2] + (b + a) / 2.;

            Gamma_pp = SQR(SQR(g)) / (64 * M_PI * SQR(mphi)) *
                       (b - a) / 2. * (nuSIaux::w_integ[0] * ((SQR(z1) - 4 * z1 + 6) / (SQR(z1) * (z1 - 2)) * log(SQR((sqrt(z1 * (z1 - 4)) + z1 - 2) / (sqrt(z1 * (z1 - 4)) - z1 + 2))) - 6 * sqrt(z1 * (z1 - 4)) / SQR(z1)) + nuSIaux::w_integ[1] * ((SQR(z2) - 4 * z2 + 6) / (SQR(z2) * (z2 - 2)) * log(SQR((sqrt(z2 * (z2 - 4)) + z2 - 2) / (sqrt(z2 * (z2 - 4)) - z2 + 2))) - 6 * sqrt(z2 * (z2 - 4)) / SQR(z2)) + nuSIaux::w_integ[2] * ((SQR(z3) - 4 * z3 + 6) / (SQR(z3) * (z3 - 2)) * log(SQR((sqrt(z3 * (z3 - 4)) + z3 - 2) / (sqrt(z3 * (z3 - 4)) - z3 + 2))) - 6 * sqrt(z3 * (z3 - 4)) / SQR(z3)));
          }

          Gamma_pp *= std::norm(U[flav][j]);
          if (majorana) // For Majorana fermions, we can scatter off neutrinos and off antineutrinos
            Gamma_pp *= 2;
        }
        tot += SQR(mphi) / (2 * mn[j]) * Gamma_pp;

        if (Gamma_s < 0 || Gamma_t_u < 0 || Gamma_tu < 0 || (Gamma_s + Gamma_t_u + Gamma_st + Gamma_su) < 0)
        {
          std::cerr << "Negative cross section when computing Gamma; for sminus/mphi^2 = "
                    << sminus << ", splus/mphi^2 = " << splus << ". The values of Gamma are as follows:" << std::endl;
          std::cerr << "Gamma_s = " << Gamma_s << std::endl;
          std::cerr << "Gamma_t_u = " << Gamma_t_u << std::endl;
          std::cerr << "Gamma_tu = " << Gamma_tu << std::endl;
          std::cerr << "Gamma_s + Gamma_t_u + Gamma_st + Gamma_su = " << Gamma_s + Gamma_t_u + Gamma_st + Gamma_su << std::endl;
          std::cerr << "Possible roundoff errors for g=" << g << ", mphi=" << mphi << ", mntot=" << mntot << std::endl;
        }
      }

      return tot;
    }

    double alphaTilde(double Em, double Ep)
    {
      /**
       * Returns
       *  \sum_{kl} \int_Em^Ep dE \int_E^Ep dE_tilde [ dsigma(E_tilde, E)/dE
       *                                             + dsigma(E_tilde, E_tilde-E)/dE]
       *                                             / |U_{flav i}|^2 |U_{flav j}|^2
       * where
       *  sigma is the regeneration cross section for the process j k -> i l
       *
       * If the final state particles are identical, there is double counting even for i != l,
       * because for final states C & D we have
       *   D=i [E_D=E], \sum_l C=l [E_C=Etilde-E]
       * + C=i [E_C=E], \sum_l D=l [E_D=Etilde-E]
       * i.e., double counting for C=D!
       */

      double Ga = scalar_width(); // Scalar decay width

      double tot = 0;
      for (int k = 0; k < 3; ++k)
      {
        /* Integration limits in terms of t/mphi^2 */
        double tplus = -2 * mn[k] * Ep / SQR(mphi);
        double tminus = -2 * mn[k] * Em / SQR(mphi);
        // If, because of a numerical coincidence, tminus==-1, shift it to avoid dividing exactly by zero
        if (fabs(tminus + 1) < 1e-7)
          tminus += tminus * 1e-6;
        // If, because of a numerical coincidence, tminus==-1, shift it to avoid dividing exactly by zero
        if (fabs(tplus + 1) < 1e-7)
          tplus += tplus * 1e-6;

        /* s-channel */
        double alphaTilde_s;
        if (fabs(tplus) < 1e-5) // We Taylor-expand atandiff to avoid roundoff errors
          alphaTilde_s = SQR(SQR(g)) / (16 * M_PI * Ga * SQR(SQR(mphi))) *
                         (2 * mphi * (1 + tminus) * (-((Ga / mphi * (1 + SQR(Ga / mphi) - 2 * tminus) * (-tminus + tplus)) / SQR(1 + SQR(Ga / mphi))) + (Ga / mphi * SQR(-tminus + tplus)) / SQR(1 + SQR(Ga / mphi))) + Ga * (log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * tplus * (tplus + 2)) -
                                                                                                                                                                                                                              log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * tminus * (tminus + 2))));
        else
          alphaTilde_s = SQR(SQR(g)) / (16 * M_PI * Ga * SQR(SQR(mphi))) *
                         (2 * mphi * (1 + tminus) * nuSIaux::atandiff(mphi * (1 + tminus) / Ga, mphi * (1 + tplus) / Ga) + Ga * (log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * tplus * (tplus + 2)) -
                                                                                                                                 log1p(SQR(mphi) / (SQR(mphi) + SQR(Ga)) * tminus * (tminus + 2))));
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alphaTilde_s *= std::norm(U[flav][k]);
        if (!majorana)
          alphaTilde_s /= 2.; // For Dirac, one of the final neutrinos is not observable
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_s;

        if (!non_resonant)
          continue;

        /* t-channel */
        double alphaTilde_t;
        if (majorana)
        {
          alphaTilde_t = SQR(SQR(g)) *
                         (1 / (16 * SQR(SQR(mphi)) * M_PI * (-1 + tminus) * tplus) *
                              ((-2 + tminus) * (tminus - tplus) -
                               (-1 + tminus) * (-2 + tplus) * (log1p(-tminus) - log1p(-tplus))) +
                          1 / (16 * SQR(SQR(mphi)) * M_PI * SQR(1 + tminus) * tplus) *
                              ((1 + tminus) * (2 + tminus) * (tminus - tplus) + (-2 * SQR(1 + tminus) + tplus + 2 * tminus * tplus) * log1p(tminus - tplus) - SQR(tminus) * tplus * log(tminus / tplus)));
          if (alphaTilde_t < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x[3], b_x[3];
            // Nodes at which the integrand will be evaluated
            double y[3], x[3][3], F[3][3];
            alphaTilde_t = 0;
            for (int i = 0; i < 3; ++i)
            {
              y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
              a_x[i] = -y[i];
              b_x[i] = -tplus;
              for (int j = 0; j < 3; ++j)
              {
                x[i][j] = (b_x[i] - a_x[i]) / 2. * nuSIaux::x_integ[j] + (b_x[i] + a_x[i]) / 2.;
                F[i][j] = SQR(y[i] / x[i][j]) / SQR(y[i] - 1) +
                          SQR((-x[i][j] - y[i]) / x[i][j]) / SQR((-x[i][j] - y[i]) - 1);
                alphaTilde_t += 1. / 4. * (b_y - a_y) * (b_x[i] - a_x[i]) * nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            }
            alphaTilde_t *= SQR(SQR(g)) / (16 * M_PI * SQR(SQR(mphi)));
          }
        }
        else
        {
          alphaTilde_t = 3. / 2. * SQR(SQR(g)) / (32 * SQR(SQR(mphi)) * M_PI * (-1 + tminus) * tplus) *
                         ((-2 + tminus) * (tminus - tplus) -
                          (-1 + tminus) * (-2 + tplus) * (log1p(-tminus) - log1p(-tplus)));
          if (alphaTilde_t < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x[3], b_x[3];
            // Nodes at which the integrand will be evaluated
            double y[3], x[3][3], F[3][3];
            alphaTilde_t = 0;
            for (int i = 0; i < 3; ++i)
            {
              y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
              a_x[i] = -y[i];
              b_x[i] = -tplus;
              for (int j = 0; j < 3; ++j)
              {
                x[i][j] = (b_x[i] - a_x[i]) / 2. * nuSIaux::x_integ[j] + (b_x[i] + a_x[i]) / 2.;
                F[i][j] = SQR(y[i] / x[i][j]) / SQR(y[i] - 1);
                alphaTilde_t += 1. / 4. * (b_y - a_y) * (b_x[i] - a_x[i]) * nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            }
            alphaTilde_t *= 3. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)));
          }
        }
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alphaTilde_t *= std::norm(U[flav][k]);
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_t;

        /* u-channel */
        double alphaTilde_u;
        if (majorana)
          alphaTilde_u = alphaTilde_t;
        else
        {
          alphaTilde_u = 1. / 2. * SQR(SQR(g)) / (32 * SQR(SQR(mphi)) * M_PI * (-1 + tminus) * tplus) *
                         ((-2 + tminus) * (tminus - tplus) -
                          (-1 + tminus) * (-2 + tplus) * (log1p(-tminus) - log1p(-tplus)));
          if (alphaTilde_u < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x[3], b_x[3];
            // Nodes at which the integrand will be evaluated
            double y[3], x[3][3], F[3][3];
            alphaTilde_u = 0;
            for (int i = 0; i < 3; ++i)
            {
              y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
              a_x[i] = -y[i];
              b_x[i] = -tplus;
              for (int j = 0; j < 3; ++j)
              {
                x[i][j] = (b_x[i] - a_x[i]) / 2. * nuSIaux::x_integ[j] + (b_x[i] + a_x[i]) / 2.;
                F[i][j] = SQR(y[i] / x[i][j]) / SQR(y[i] - 1);
                alphaTilde_u += 1. / 4. * (b_y - a_y) * (b_x[i] - a_x[i]) * nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            }
            alphaTilde_u *= 1. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)));
          }
          // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
          alphaTilde_u *= std::norm(U[flav][k]);
        }
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_u;

        /* t-u interference */
        double alphaTilde_tu;
        if (majorana)
        {
          double dilog_combi;
          if (-tplus < 1e-2 && -tminus < 1e-2)
          {
            double delta = tplus / tminus;
            dilog_combi = -(((-1 + delta) * tplus * log(-2 * tplus)) / delta) - ((-1 + delta) * SQR(tplus) * (-2 + delta + delta * log(2) + log(-2 / tplus) - delta * log(-tplus))) / (2. * SQR(delta)) +
                          (CUB(tplus) * (8 - 30 * delta + 21 * SQR(delta) + CUB(delta) - 8 * CUB(delta) * log(2) + log(256) + 8 * log(-tplus) - 8 * CUB(delta) * log(-tplus))) /
                              (24. * CUB(delta)) +
                          (SQR(SQR(tplus)) * (-32 + 56 * delta - 51 * SQR(delta) + 30 * CUB(delta) - 3 * SQR(SQR(delta)) + log(4096) - SQR(SQR(delta)) * log(4096) -
                                              12 * log(-tplus) + 12 * SQR(SQR(delta)) * log(-tplus))) /
                              (48. * SQR(SQR(delta)));
          }
          else if (-tplus > 1e2 && -tminus > 1e2)
          {
            double delta = tplus / tminus;
            dilog_combi = (-2 * (-1 + delta) * log((-1 + delta) / delta)) / tplus - (2 * (-1 + log(-(delta / ((-1 + delta) * tplus))))) / SQR(tplus) +
                          (-6 + 4 * delta + SQR(delta) - 2 * CUB(delta) - 8 * log((-1 + delta) / delta) + 8 * delta * log((-1 + delta) / delta) + 2 * CUB(delta) * log((-1 + delta) / delta) -
                           2 * SQR(SQR(delta)) * log((-1 + delta) / delta) - 6 * log(-tplus) + 6 * delta * log(-tplus)) /
                              (3. * (-1 + delta) * CUB(tplus)) +
                          (8 - 12 * delta + 3 * SQR(delta) + 12 * log((-1 + delta) / delta) - 24 * delta * log((-1 + delta) / delta) + 12 * SQR(delta) * log((-1 + delta) / delta) + 12 * log(-tplus) -
                           24 * delta * log(-tplus) + 12 * SQR(delta) * log(-tplus)) /
                              (3. * SQR(-1 + delta) * SQR(SQR(tplus)));
          }
          else
            dilog_combi = gsl_sf_dilog(1 + 1 / (-2 + tplus)) - gsl_sf_dilog((-1 + tminus) / (-2 + tplus)) + gsl_sf_dilog(1 + (1 + tminus - tplus) / tplus) - gsl_sf_dilog(1 + 1 / tplus);

          alphaTilde_tu = SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)) * (1 + tminus) * tplus) *
                          (2 * (2 * (1 + tminus) * (tminus - tplus) - 2 * (1 + tminus) * tplus * atanh(1 / (1 - tplus)) * atanh((tminus - tplus) / (-2 + tminus + tplus)) +
                                tminus * tplus * (-log1p(-tminus) + log1p(-tplus)) + (1 + tminus) * (log1p(-tminus) - log1p(-tplus) - log1p(tminus - tplus)) +
                                tplus * (-log1p(-tminus) + log1p(-tplus) + log1p(tminus - tplus)) - tminus * tplus * log(tminus / tplus)) +
                           (1 + tminus) * tplus * ((-SQR(log1p(-tminus)) + SQR(log1p(-tplus))) / 2. + nuSIaux::dilog1over1mdiff(tplus, tminus)) -
                           (1 + tminus) * tplus * (nuSIaux::dilog1pdiff(tminus, tplus) + dilog_combi));

          if (alphaTilde_tu < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x[3], b_x[3];
            // Nodes at which the integrand will be evaluated
            double y[3], x[3][3], F[3][3];
            alphaTilde_tu = 0;
            for (int i = 0; i < 3; ++i)
            {
              y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
              a_x[i] = -y[i];
              b_x[i] = -tplus;
              for (int j = 0; j < 3; ++j)
              {
                x[i][j] = (b_x[i] - a_x[i]) / 2. * nuSIaux::x_integ[j] + (b_x[i] + a_x[i]) / 2.;
                F[i][j] = 2 * y[i] * (-y[i] - x[i][j]) / SQR(x[i][j]) / ((y[i] - 1) * (-y[i] - x[i][j] - 1));
                alphaTilde_tu += 1. / 4. * (b_y - a_y) * (b_x[i] - a_x[i]) * nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            }
            alphaTilde_tu *= SQR(SQR(g)) / (16 * M_PI * SQR(SQR(mphi)));
          }
        }
        else
          alphaTilde_tu = 0;
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alphaTilde_tu *= std::norm(U[flav][k]);
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_tu;

        /* s-t interference */
        double ga_red = Ga / mphi; // Reduced decay width
        // Auxiliary complex numbers
        double _Complex z1 = (-I * (-1 + tminus)) / (2 * I + ga_red);
        double z2 = 1 / (1 + tminus);
        double _Complex z3 = 1 / (2 - I * ga_red + tminus);
        double _Complex z4 = (1 + tminus - tplus) / (2 - I * ga_red + tminus);
        double _Complex z5 = (-I * (-1 + tplus)) / (2 * I + ga_red);
        double z6 = 1 - tplus / (1 + tminus);
        double z7 = 1 - tminus;
        double z8 = 1 - tplus;
        // Compute dilogarithms
        double _Complex dilogdiff_z7z8;
        double _Complex dilogdiff_z5z1;
        double _Complex dilogdiff_z2z6;
        double _Complex dilogdiff_z4z3;

        if (-tplus < 1e-5)
        { // We Taylor-expand dilogdiff to be fast and avoid roundoff errors
          double delta = tplus / tminus;
          dilogdiff_z7z8 = tminus * (-1 + clog(tminus)) + (SQR(tminus) * (-1 + 2 * clog(tminus))) / 4. -
                           (tplus * (-1 + clog(tplus)) + (SQR(tplus) * (-1 + 2 * clog(tplus))) / 4.);
          dilogdiff_z5z1 = (-tminus + tplus) * clog(1 - I / (2 * I + ga_red)) +
                           ((-SQR(tminus) + SQR(tplus)) * (I * (1 + clog(1 - I / (2 * I + ga_red))) +
                                                           clog(1 - I / (2 * I + ga_red)) * ga_red)) /
                               (2. * (I + ga_red));
          dilogdiff_z2z6 = (tplus * (-1 + delta - clog(delta) + clog(tplus) - delta * clog(tplus))) / delta +
                           (SQR(tplus) * (-1 + SQR(delta) + 2 * clog(delta) - 2 * clog(tplus) + 4 * delta * clog(tplus) - 2 * SQR(delta) * clog(tplus))) / (4. * SQR(delta)) +
                           (CUB(tplus) * (7 - 9 * delta + 2 * CUB(delta) - 6 * clog(delta) + 6 * clog(tplus) - 18 * delta * clog(tplus) + 18 * SQR(delta) * clog(tplus) -
                                          6 * CUB(delta) * clog(tplus))) /
                               (18. * CUB(delta));
          dilogdiff_z4z3 = ((-1 + delta) * tplus * clog((I + ga_red) / (2 * I + ga_red))) / delta +
                           ((-1 + delta) * SQR(tplus) * (I * ((1 + delta) / (I + ga_red) - 2 / (2 * I + ga_red)) + (-1 + delta) * clog((I + ga_red) / (2 * I + ga_red)))) / (2. * SQR(delta));
        }
        else
        {
          dilogdiff_z7z8 = nuSIaux::dilogdiff_complex(z7, z8);
          dilogdiff_z5z1 = nuSIaux::dilogdiff_complex(z5, z1);
          dilogdiff_z2z6 = nuSIaux::dilogdiff_complex(z2, z6);
          dilogdiff_z4z3 = nuSIaux::dilogdiff_complex(z4, z3);
        }

        double alphaTilde_st;
        if (majorana)
          alphaTilde_st = SQR(SQR(g)) / (32 * M_PI * (1 + SQR(ga_red)) * SQR(SQR(mphi))) *
                          (2 * M_PI * carg(-1 + I * ga_red - tminus) - 2 * M_PI * carg(-1 + I * ga_red - tplus) + 2 * ga_red * (cimag(dilogdiff_z5z1) + cimag(dilogdiff_z2z6) + cimag(dilogdiff_z4z3)) - 2 * (creal(dilogdiff_z5z1) + creal(dilogdiff_z2z6) + creal(dilogdiff_z4z3) + creal(dilogdiff_z7z8)) - carg((ga_red + I * (1 + tminus)) / (2 * I + ga_red)) * (2 * M_PI + 2 * ga_red * log1p(-tminus)) + carg((ga_red + I * (1 + tplus)) / (2 * I + ga_red)) * (2 * M_PI + 2 * ga_red * log1p(-tplus)) + (carg(-1 + I * ga_red - tminus) - carg(-1 + I * ga_red - tplus)) * (4 * ga_red * tminus + 2 * ga_red * log1p(-tminus)) + 2 * ga_red * (carg(1 + tminus) - carg(2 - I * ga_red + tminus) + carg(1 - I * ga_red + tplus)) * log1p(tminus - tplus) + log(4 + SQR(ga_red)) * (log1p(-tplus) - log1p(-tminus)) + log(SQR(ga_red) + SQR(2 + tminus)) * log1p(tminus - tplus) - 2 * log1p(-tminus) * log(-tplus) - 2 * ga_red * M_PI * (log(SQR(tplus)) + log1p(tminus - tplus)) + 2 * ga_red * M_PI * log(SQR(tplus)) + 4 * tminus * log(tminus / tplus) + (-log1p(-tplus) + log1p(-tminus) - log1p(tminus - tplus)) * (log1p(SQR(1 + tplus) / SQR(ga_red)) + 2 * log(ga_red)) - log1p(tminus - tplus) * log1p(SQR(tminus) + 2 * tminus) + 2 * (SQR(ga_red) + tminus) * (log1p(SQR(1 + tplus) / SQR(ga_red)) - log1p(SQR(1 + tminus) / SQR(ga_red))) + 2 * (log(-tplus) * (log1p(-tplus) + log1p(tminus - tplus)) + (log1p(SQR(1 + tplus) / SQR(ga_red)) - log1p(SQR(1 + tminus) / SQR(ga_red)))));
        else
          alphaTilde_st = SQR(SQR(g)) / (32 * M_PI * (1 + SQR(ga_red)) * SQR(SQR(mphi))) *
                          (ga_red * cimag(dilogdiff_z5z1) - 2 * (creal(dilogdiff_z5z1 + dilogdiff_z7z8)) + 2 * carg((ga_red + I * (1 + tminus)) / (2 * I + ga_red)) * (-M_PI - ga_red * log1p(-tminus)) + 2 * carg(-1 + I * ga_red - tminus) * (M_PI + ga_red * tminus + ga_red * log1p(-tminus)) - 2 * carg(-1 + I * ga_red - tplus) * (M_PI + ga_red * tminus + ga_red * log1p(-tminus)) + 2 * carg((ga_red + I * (1 + tplus)) / (2 * I + ga_red)) * (M_PI + ga_red * log1p(-tplus)) - 2 * log1p(-tminus) * log(-tplus) + 2 * tminus * log(tminus / tplus) + 2 * log1p(-tplus) * log(-tplus) + (log1p(-tplus) - log1p(-tminus)) * (log(4 + SQR(ga_red)) - 2 * log(ga_red) - log1p(SQR(1 + tplus) / SQR(ga_red))) + (1 + tminus + SQR(ga_red)) * (log1p(SQR(1 + tplus) / SQR(ga_red)) - log1p(SQR(1 + tminus) / SQR(ga_red))));

        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alphaTilde_st *= std::norm(U[flav][k]);
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_st;

        /* s-u interference (only for Majorana fermions)*/
        double alphaTilde_su = 0;
        if (majorana)
          alphaTilde_su = alphaTilde_st;
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_su;

        /* Double scalar production */
        double alphaTilde_pp = 0;
        if ((-tplus > 4) && phiphi)
        {
          if (-tplus < 1e4)
            alphaTilde_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * spl_alphaTilde_phiphi.f_eval({-tplus, log10(tplus / tminus)});
          else
          { // Use a Taylor expansion for the non interpolated region
            alphaTilde_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * (6 * tminus * log(-tminus) - tplus * SQR(log(-tminus)) + 2 * (-8 * tminus + 8 * tplus + 4 * tplus * log(-tminus) + log(tminus - tplus) * (tminus - tplus - tplus * log(tminus / tplus))) - 2 * (2 * tminus + 5 * tplus) * log(-tplus) + tplus * SQR(log(-tplus)) - 2 * tplus * gsl_sf_dilog(1 - tminus / tplus)) / (128. * M_PI * tplus);
          }

          alphaTilde_pp *= std::norm(U[flav][k]);
          if (majorana) // For Majorana fermions, we can scatter off neutrinos and off antineutrinos
            alphaTilde_pp *= 2;

          alphaTilde_pp *= 2; // Each scattering generates at least 2 neutrinos
          if (majorana)       // For Majorana fermions, all final states are observable
            alphaTilde_pp *= 2;
        }
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alphaTilde_pp;

        if (alphaTilde_s < 0 || alphaTilde_t < 0 || alphaTilde_u < 0 || alphaTilde_tu / SQR(SQR(g / mphi)) < -1e-11
            /* In the s-t interference, there can be roundoff errors from computing dilog(1+x+delta) - dilog(1+x) for delta<<<
               To remove them, a dedicated Taylor expansion would be necessary. Thus, we will ignore them as long as they are numerically irrelevant
               I.e., we'll impose that the cross section is positive **within a reasonable numerical precision**
            */
            || (alphaTilde_st + alphaTilde_t + alphaTilde_s) / SQR(SQR(g / mphi)) < -1e-11 || (alphaTilde_su + alphaTilde_u + alphaTilde_s) / SQR(SQR(g / mphi)) < -1e-11)
        {
          std::cerr << "Negative cross section when computing alphaTilde; for tminus/mphi^2 = "
                    << tminus << ", tplus/mphi^2 = " << tplus << ". The values of alphaTilde are as follows:" << std::endl;
          std::cerr << "alphaTilde_s = " << alphaTilde_s << std::endl;
          std::cerr << "alphaTilde_t = " << alphaTilde_t << std::endl;
          std::cerr << "alphaTilde_u = " << alphaTilde_u << std::endl;
          std::cerr << "alphaTilde_tu = " << alphaTilde_tu << std::endl;
          std::cerr << "alphaTilde_st + alphaTilde_s + alphaTilde_t = " << alphaTilde_st + alphaTilde_s + alphaTilde_t << std::endl;
          std::cerr << "alphaTilde_su + alphaTilde_s + alphaTilde_u = " << alphaTilde_su + alphaTilde_s + alphaTilde_u << std::endl;
          std::cerr << "Possible roundoff errors for g=" << g << ", mphi=" << mphi << ", mntot=" << mntot << std::endl;
        }
      }

      return tot;
    }

    double alpha(double Em, double Ep, double Em_prime, double Ep_prime)
    {
      /**
       * Returns
       *  \sum_{kl} \int_Em^Ep dE \int_Em_prime^Ep_prime dE_tilde [ dsigma(E_tilde, E)/dE
       *                                                          + dsigma(E_tilde, (E_tilde-E))/dE] / |U_{flav i}|^2 |U_{flav j}|^2
       * where
       *  sigma is the regeneration cross section for the process j k -> i l
       */

      double Ga = scalar_width(); // Scalar decay width

      double tot = 0;
      for (int k = 0; k < 3; ++k)
      {
        /* Integration limits in terms of t/mphi^2 and s/mphi^2 */
        double tplus = -2 * mn[k] * Ep / SQR(mphi);
        double tminus = -2 * mn[k] * Em / SQR(mphi);
        double splus_prime = 2 * mn[k] * Ep_prime / SQR(mphi);
        double sminus_prime = 2 * mn[k] * Em_prime / SQR(mphi);
        // If, because of a numerical coincidence, tminus==-1, shift it to avoid dividing exactly by zero
        if (fabs(tminus + 1) < 1e-7)
          tminus += tminus * 1e-6;
        // If, because of a numerical coincidence, tminus==-1, shift it to avoid dividing exactly by zero
        if (fabs(tplus + 1) < 1e-7)
          tplus += tplus * 1e-6;

        /* s-channel */
        double alpha_s;
        if (splus_prime < 1e-5) // We Taylor-expand atandiff to avoid roundoff errors
          alpha_s = SQR(SQR(g)) / (8 * M_PI * Ga * CUB(mphi)) * (tminus - tplus) * ((Ga / mphi * (1 + SQR(Ga / mphi) + 2 * sminus_prime)) / SQR(1 + SQR(Ga / mphi)) * (splus_prime - sminus_prime) + (Ga / mphi) / SQR(1 + SQR(Ga / mphi)) * SQR(splus_prime - sminus_prime));
        else
          alpha_s = SQR(SQR(g)) / (8 * M_PI * Ga * CUB(mphi)) * (tminus - tplus) * nuSIaux::atandiff(mphi * (splus_prime - 1) / Ga, mphi * (sminus_prime - 1) / Ga);
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alpha_s *= std::norm(U[flav][k]);
        if (!majorana)
          alpha_s /= 2.; // For Dirac, one of the final neutrinos is not observable

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_s;

        if (!non_resonant)
          continue;

        /* t-channel */
        double alpha_t;
        if (majorana)
        {
          alpha_t = SQR(SQR(g)) / (sminus_prime * splus_prime * 16 * M_PI * SQR(SQR(mphi))) * (-((sminus_prime - splus_prime) * (3 + 2 * tminus * (-1 + tplus) - 2 * tplus) * (tminus - tplus)) / ((-1 + tminus) * (-1 + tplus)) + 2 * (sminus_prime * splus_prime * (-tminus + tplus) * log(sminus_prime) + sminus_prime * splus_prime * (tminus - tplus) * log(splus_prime) - sminus_prime * splus_prime * log1p(sminus_prime + tminus) - sminus_prime * splus_prime * tplus * log1p(sminus_prime + tminus) + sminus_prime * splus_prime * log1p(splus_prime + tminus) + sminus_prime * splus_prime * tplus * log1p(splus_prime + tminus) - splus_prime * log(((1 + sminus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + sminus_prime + tplus))) - splus_prime * tminus * log(((1 + sminus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + sminus_prime + tplus))) - splus_prime * tplus * log(((1 + sminus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + sminus_prime + tplus))) - splus_prime * tminus * tplus * log(((1 + sminus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + sminus_prime + tplus))) + sminus_prime * splus_prime * log(1 + sminus_prime + tplus) + sminus_prime * splus_prime * tminus * log1p(sminus_prime + tplus) + sminus_prime * log(((1 + splus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + splus_prime + tplus))) + sminus_prime * tminus * log(((1 + splus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + splus_prime + tplus))) + sminus_prime * tplus * log(((1 + splus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + splus_prime + tplus))) + sminus_prime * tminus * tplus * log(((1 + splus_prime + tminus) * (-1 + tplus)) / ((-1 + tminus) * (1 + splus_prime + tplus))) - sminus_prime * splus_prime * log(1 + splus_prime + tplus) - sminus_prime * splus_prime * tminus * log1p(splus_prime + tplus)) / ((1 + tminus) * (1 + tplus)) - ((sminus_prime * splus_prime * log((sminus_prime * (1 + splus_prime + tminus)) / (splus_prime * (1 + sminus_prime + tminus)))) / SQR(1 + tminus) + (((sminus_prime - splus_prime) * (tminus - tplus) * (1 + tplus)) / (1 + tminus) - sminus_prime * splus_prime * log((sminus_prime * (1 + splus_prime + tplus)) / (splus_prime * (1 + sminus_prime + tplus)))) / SQR(1 + tplus)));

          if (alpha_t < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x = sminus_prime, b_x = splus_prime;
            // Nodes at which the integrand will be evaluated
            double y[3], x[3], F[3][3];
            alpha_t = 0;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j)
              {
                y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
                x[j] = (b_x - a_x) / 2. * nuSIaux::x_integ[j] + (b_x + a_x) / 2.;
                F[i][j] = SQR(y[i] / x[j]) / SQR(y[i] - 1) +
                          SQR((-x[j] - y[i]) / x[j]) / SQR((-x[j] - y[i]) - 1);
                alpha_t += nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            alpha_t *= 1. / 4. * (b_y - a_y) * (b_x - a_x);

            alpha_t *= SQR(SQR(g)) / (16 * M_PI * SQR(SQR(mphi)));
          }
        }
        else
        {
          alpha_t = 3. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)) * sminus_prime * splus_prime * (-1 + tminus) * (-1 + tplus)) *
                    (sminus_prime - splus_prime) * (-((tminus - tplus) * (2 + tminus * (-1 + tplus) - tplus)) - 2 * (-1 + tminus) * (-1 + tplus) * (log1p(-tminus) - log1p(-tplus)));

          if (alpha_t < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x = sminus_prime, b_x = splus_prime;
            // Nodes at which the integrand will be evaluated
            double y[3], x[3], F[3][3];
            alpha_t = 0;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j)
              {
                y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
                x[j] = (b_x - a_x) / 2. * nuSIaux::x_integ[j] + (b_x + a_x) / 2.;
                F[i][j] = SQR(y[i] / x[j]) / SQR(y[i] - 1);
                alpha_t += nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            alpha_t *= 1. / 4. * (b_y - a_y) * (b_x - a_x);

            alpha_t *= 3. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)));
          }
        }
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alpha_t *= std::norm(U[flav][k]);

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_t;

        /* u-channel */
        double alpha_u;
        if (majorana)
          alpha_u = alpha_t;
        else
        {
          alpha_u = 1. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)) * sminus_prime * splus_prime * (-1 + tminus) * (-1 + tplus)) *
                    (sminus_prime - splus_prime) * (-((tminus - tplus) * (2 + tminus * (-1 + tplus) - tplus)) - 2 * (-1 + tminus) * (-1 + tplus) * (log1p(-tminus) - log1p(-tplus)));

          if (alpha_u < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x = sminus_prime, b_x = splus_prime;
            // Nodes at which the integrand will be evaluated
            double y[3], x[3], F[3][3];
            alpha_u = 0;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j)
              {
                y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
                x[j] = (b_x - a_x) / 2. * nuSIaux::x_integ[j] + (b_x + a_x) / 2.;
                F[i][j] = SQR(y[i] / x[j]) / SQR(y[i] - 1);
                alpha_u += nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            alpha_u *= 1. / 4. * (b_y - a_y) * (b_x - a_x);

            alpha_u *= 1. / 2. * SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)));
          }

          // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
          alpha_u *= std::norm(U[flav][k]);
        }

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_u;

        /* t-u interference */
        double alpha_tu;
        if (majorana)
        {
          double FCTR_tplus;
          if (tplus < -1)
            FCTR_tplus = gsl_sf_dilog((1 + sminus_prime + tplus) / sminus_prime) - gsl_sf_dilog((1 + splus_prime + tplus) / splus_prime);
          else
            FCTR_tplus = -gsl_sf_dilog(sminus_prime / (1 + sminus_prime + tplus)) + gsl_sf_dilog(splus_prime / (1 + splus_prime + tplus)) -
                         0.5 * (SQR(log((1 + sminus_prime + tplus) / sminus_prime)) - SQR(log((1 + splus_prime + tplus) / splus_prime)));
          double FCTR_tminus;
          if (tminus < -1)
            FCTR_tminus = -gsl_sf_dilog((1 + sminus_prime + tminus) / sminus_prime) + gsl_sf_dilog((1 + splus_prime + tminus) / splus_prime);
          else
            FCTR_tminus = gsl_sf_dilog(sminus_prime / (1 + sminus_prime + tminus)) - gsl_sf_dilog(splus_prime / (1 + splus_prime + tminus)) +
                          0.5 * (SQR(log((1 + sminus_prime + tminus) / sminus_prime)) - SQR(log((1 + splus_prime + tminus) / splus_prime)));

          double log1p_abs_tplus = (tplus > -1) ? log1p(tplus) : log(-1 - tplus);     // log|1+tplus|
          double log1p_abs_tminus = (tminus > -1) ? log1p(tminus) : log(-1 - tminus); // log|1+tminus|

          alpha_tu = SQR(SQR(g)) / (32 * M_PI * SQR(SQR(mphi)) * sminus_prime * splus_prime * (1 + tminus) * (1 + tplus)) *
                     (-4 * (sminus_prime - splus_prime) * (1 + tminus) * (tminus - tplus) * (1 + tplus) +
                      2 * sminus_prime * splus_prime * tplus * (log(sminus_prime / splus_prime) - log1p(sminus_prime + tminus) + log1p(splus_prime + tminus)) +
                      2 * splus_prime * (1 + tminus) * (1 + tplus) * (log1p(-tminus) - log1p(sminus_prime + tminus) - log1p(-tplus) + log1p(sminus_prime + tplus)) -
                      2 * sminus_prime * (1 + tminus) * (1 + tplus) * (log1p(-tminus) - log1p(splus_prime + tminus) - log1p(-tplus) + log1p(splus_prime + tplus)) +
                      2 * sminus_prime * splus_prime * (-log1p(sminus_prime + tminus) + log1p(splus_prime + tminus) + log1p(sminus_prime + tplus) - log1p(splus_prime + tplus)) +
                      sminus_prime * splus_prime * (1 + tminus) * (1 + tplus) * (log((2 + sminus_prime) / sminus_prime) * (log(splus_prime) + log1p(sminus_prime + tplus)) - log((2 + splus_prime) / splus_prime) * (log(sminus_prime) + log1p(splus_prime + tplus)) + log1p(-tplus) * (log(sminus_prime / splus_prime) - log1p(sminus_prime + tplus) + log1p(splus_prime + tplus))) +
                      sminus_prime * splus_prime * (1 + tminus) * (1 + tplus) * ((log(splus_prime) + log1p(sminus_prime + tminus)) * (log(sminus_prime / (2 + sminus_prime)) + log1p(-tminus) - log1p_abs_tminus) + (log(sminus_prime) + log1p(splus_prime + tminus)) * (log((2 + splus_prime) / splus_prime) - log1p(-tminus) + log1p_abs_tminus)) +
                      sminus_prime * splus_prime * (log(splus_prime / sminus_prime) + log1p(sminus_prime + tplus) - log1p(splus_prime + tplus)) * (2 * tminus + (1 + tminus) * (1 + tplus) * log1p_abs_tplus) +
                      sminus_prime * splus_prime * (1 + tminus) * (1 + tplus) * (gsl_sf_dilog((1 + sminus_prime + tminus) / (2 + sminus_prime)) - gsl_sf_dilog((1 + splus_prime + tminus) / (2 + splus_prime)) - gsl_sf_dilog((1 + sminus_prime + tplus) / (2 + sminus_prime)) + gsl_sf_dilog((1 + splus_prime + tplus) / (2 + splus_prime))) +
                      sminus_prime * splus_prime * (1 + tminus) * (1 + tplus) * (FCTR_tplus + FCTR_tminus));

          if (alpha_tu < 0)
          { // Roundoff errors! Compute the integral numerically
            double a_y = tplus, b_y = tminus, a_x = sminus_prime, b_x = splus_prime;
            // Nodes at which the integrand will be evaluated
            double y[3], x[3], F[3][3];
            double alpha_tu = 0;
            for (int i = 0; i < 3; ++i)
              for (int j = 0; j < 3; ++j)
              {
                y[i] = (b_y - a_y) / 2. * nuSIaux::x_integ[i] + (b_y + a_y) / 2.;
                x[j] = (b_x - a_x) / 2. * nuSIaux::x_integ[j] + (b_x + a_x) / 2.;
                F[i][j] = 2 * y[i] * (-y[i] - x[j]) / SQR(x[j]) / ((y[i] - 1) * (-y[i] - x[j] - 1));
                alpha_tu += nuSIaux::w_integ[i] * nuSIaux::w_integ[j] * F[i][j];
              }
            alpha_tu *= 1. / 4. * (b_y - a_y) * (b_x - a_x);

            alpha_tu *= SQR(SQR(g)) / (16 * M_PI * SQR(SQR(mphi)));
          }
        }
        else
          alpha_tu = 0.;
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alpha_tu *= std::norm(U[flav][k]);

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_tu;

        /* s-t interference */
        double alpha_st;
        double ga_red = Ga / mphi; // Reduced decay width
        // Auxiliary complex numbers
        double z1 = (1 + sminus_prime + tminus) / (1 + tminus);
        double _Complex z2 = (1 + sminus_prime + tminus) / (2 - I * ga_red + tminus);
        double z3 = (1 + splus_prime + tminus) / (1 + tminus);
        double _Complex z4 = (1 + splus_prime + tminus) / (2 - I * ga_red + tminus);
        double z5 = (1 + sminus_prime + tplus) / (1 + tplus);
        double _Complex z6 = (1 + sminus_prime + tplus) / (2 - I * ga_red + tplus);
        double z7 = (1 + splus_prime + tplus) / (1 + tplus);
        double _Complex z8 = (1 + splus_prime + tplus) / (2 - I * ga_red + tplus);
        // Compute dilogarithms
        gsl_sf_result dilog_z1_re, dilog_z1_im, dilog_z2_re, dilog_z2_im,
            dilog_z3_re, dilog_z3_im, dilog_z4_re, dilog_z4_im,
            dilog_z5_re, dilog_z5_im, dilog_z6_re, dilog_z6_im,
            dilog_z7_re, dilog_z7_im, dilog_z8_re, dilog_z8_im;
        gsl_sf_complex_dilog_xy_e(z1, 0, &dilog_z1_re, &dilog_z1_im);
        gsl_sf_complex_dilog_xy_e(creal(z2), cimag(z2), &dilog_z2_re, &dilog_z2_im);
        gsl_sf_complex_dilog_xy_e(z3, 0, &dilog_z3_re, &dilog_z3_im);
        gsl_sf_complex_dilog_xy_e(creal(z4), cimag(z4), &dilog_z4_re, &dilog_z4_im);
        gsl_sf_complex_dilog_xy_e(z5, 0, &dilog_z5_re, &dilog_z5_im);
        gsl_sf_complex_dilog_xy_e(creal(z6), cimag(z6), &dilog_z6_re, &dilog_z6_im);
        gsl_sf_complex_dilog_xy_e(z7, 0, &dilog_z7_re, &dilog_z7_im);
        gsl_sf_complex_dilog_xy_e(creal(z8), cimag(z8), &dilog_z8_re, &dilog_z8_im);

        if (majorana)
        {
          alpha_st = SQR(SQR(g)) / (32 * M_PI * (1 + SQR(ga_red)) * SQR(SQR(mphi))) *
                     (2 * ga_red * (dilog_z1_im.val - dilog_z2_im.val - dilog_z3_im.val + dilog_z4_im.val - dilog_z5_im.val + dilog_z6_im.val + dilog_z7_im.val - dilog_z8_im.val) - 2 * (dilog_z1_re.val - dilog_z2_re.val - dilog_z3_re.val + dilog_z4_re.val - dilog_z5_re.val + dilog_z6_re.val + dilog_z7_re.val - dilog_z8_re.val) + 2 * ga_red * (carg(-(1 / (1 + tminus))) - carg(-((-1 + I * ga_red + sminus_prime) / (2 - I * ga_red + tminus)))) * log1p(sminus_prime + tminus) - 2 * ga_red * (carg(-(1 / (1 + tminus))) - carg(-((-1 + I * ga_red + splus_prime) / (2 - I * ga_red + tminus)))) * log1p(splus_prime + tminus) + 2 * ga_red * (carg(-(1 / (1 + tplus))) - carg(-((-1 + I * ga_red + splus_prime) / (2 - I * ga_red + tplus)))) * log1p(splus_prime + tplus) - 2 * ga_red * (carg(-(1 / (1 + tplus))) - carg(-((-1 + I * ga_red + sminus_prime) / (2 - I * ga_red + tplus)))) * log1p(sminus_prime + tplus) + 2 * (ga_red * carg(-1 + I * ga_red + sminus_prime) - ga_red * carg(-1 + I * ga_red + splus_prime) + log1p(SQR(-1 + splus_prime) / SQR(ga_red)) / 2. - log1p(SQR(-1 + sminus_prime) / SQR(ga_red)) / 2. + log(sminus_prime) - log(splus_prime)) * (2 * (tminus - tplus) + (log1p(-tminus) - log1p(-tplus))) + log1p(sminus_prime + tminus) * (log1p(SQR(-1 + sminus_prime) / SQR(ga_red)) - log1p(SQR(2 + tminus) / SQR(ga_red)) - 2 * (log(sminus_prime) - log(fabs(1 + tminus)))) - log1p(splus_prime + tminus) * (log1p(SQR(-1 + splus_prime) / SQR(ga_red)) - log1p(SQR(2 + tminus) / SQR(ga_red)) - 2 * (log(splus_prime) - log(fabs(1 + tminus)))) - log1p(sminus_prime + tplus) * (log1p(SQR(-1 + sminus_prime) / SQR(ga_red)) - log1p(SQR(2 + tplus) / SQR(ga_red)) - 2 * (log(sminus_prime) - log(fabs(1 + tplus)))) + log1p(splus_prime + tplus) * (log1p(SQR(-1 + splus_prime) / SQR(ga_red)) - log1p(SQR(2 + tplus) / SQR(ga_red)) - 2 * (log(splus_prime) - log(fabs(1 + tplus)))));
        }
        else
        {
          alpha_st = SQR(SQR(g)) / (32 * M_PI * (1 + SQR(ga_red)) * SQR(SQR(mphi))) *
                     ((2 * ga_red * carg(-1 + I * ga_red + sminus_prime) - 2 * ga_red * carg(-1 + I * ga_red + splus_prime) + 2 * log(sminus_prime) - 2 * log(splus_prime) + log1p(SQR(-1 + splus_prime) / SQR(ga_red)) - log1p(SQR(-1 + sminus_prime) / SQR(ga_red))) *
                      (tminus - tplus + log1p(-tminus) - log1p(-tplus)));
        }
        // Prefactor: |U_{flav i}|^2 |U_{flav j}|^2 |U_{flav k}|^2 \sum_{l} |U_{flav l}|^2
        alpha_st *= std::norm(U[flav][k]);

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_st;

        /* s-u interference */
        double alpha_su = 0.;
        if (majorana)
          alpha_su = alpha_st;

        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_su;

        /* Double scalar production */
        double alpha_pp = 0;
        if ((sminus_prime > 4) && phiphi)
        {
          if (sminus_prime < 1e4)
          {
            double delta = splus_prime / sminus_prime;
            alpha_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * fabs(spl_alpha_phiphi.f_eval({sminus_prime, log(-sminus_prime / tminus) / log(delta) * 1.0001, log10(delta)}));
          }
          else
          { // Use a Taylor expansion for the non interpolated region
            if (tminus < -1)
              alpha_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * ((-sminus_prime + splus_prime) * ((tminus - tplus) * (splus_prime * (-2 + tminus + tplus) + sminus_prime * (-2 - 24 * splus_prime + tminus + tplus)) + 4 * (-(splus_prime * (1 + tminus)) + sminus_prime * (-1 + 2 * splus_prime + (-1 + splus_prime) * tminus)) * log(-1 - tminus) + 2 * (3 * splus_prime + sminus_prime * (3 + 4 * splus_prime)) * tminus * log(-tminus) + 4 * (splus_prime + splus_prime * tplus + sminus_prime * (1 + tplus - splus_prime * (2 + tplus))) * log(-1 - tplus) - 2 * (3 * splus_prime + sminus_prime * (3 + 4 * splus_prime)) * tplus * log(-tplus)) + 2 * SQR(sminus_prime) * log(splus_prime) * ((3 + 2 * splus_prime) * (tminus - tplus) + 2 * SQR(splus_prime) * ((-1 - tminus) * log(-1 - tminus) + tminus * log(-tminus) + (1 + tplus) * log(-1 - tplus) - tplus * log(-tplus))) + 2 * SQR(splus_prime) * log(sminus_prime) * ((-3 - 2 * sminus_prime) * (tminus - tplus) + 2 * SQR(sminus_prime) * ((1 + tminus) * log(-1 - tminus) - tminus * log(-tminus) - (1 + tplus) * log(-1 - tplus) + tplus * log(-tplus)))) / (256. * M_PI * SQR(sminus_prime) * SQR(splus_prime));
            else if (tplus < -1)
              alpha_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * ((2 * SQR(sminus_prime) * log(splus_prime) * ((1 + tplus) * (-3 - 2 * splus_prime + 2 * SQR(splus_prime) * log(-1 - tplus)) - 2 * SQR(splus_prime) * tplus * log(-tplus)) + (sminus_prime - splus_prime) * ((1 + tplus) * (-3 * (sminus_prime + splus_prime + 8 * sminus_prime * splus_prime) + (sminus_prime + splus_prime) * tplus) + 4 * (-(splus_prime * (1 + tplus)) + sminus_prime * (-1 + 2 * splus_prime + (-1 + splus_prime) * tplus)) * log(-1 - tplus) + 2 * (3 * splus_prime + sminus_prime * (3 + 4 * splus_prime)) * tplus * log(-tplus)) + 2 * SQR(splus_prime) * log(sminus_prime) * ((3 + 2 * sminus_prime) * (1 + tplus) + 2 * SQR(sminus_prime) * (-((1 + tplus) * log(-1 - tplus)) + tplus * log(-tplus)))) / (256. * M_PI * SQR(sminus_prime) * SQR(splus_prime)) + (-1 - tminus) * (-6 * sminus_prime + 6 * splus_prime - 2 * (-2 + sminus_prime) * splus_prime * log(sminus_prime) + sminus_prime * splus_prime * SQR(log(sminus_prime)) + 2 * sminus_prime * (-2 + splus_prime) * log(splus_prime) - sminus_prime * splus_prime * SQR(log(splus_prime))) / (128. * M_PI * sminus_prime * splus_prime));
            else
              alpha_pp = SQR(SQR(g)) / SQR(SQR(mphi)) * (tplus - tminus) * (-6 * sminus_prime + 6 * splus_prime - 2 * (-2 + sminus_prime) * splus_prime * log(sminus_prime) + sminus_prime * splus_prime * SQR(log(sminus_prime)) + 2 * sminus_prime * (-2 + splus_prime) * log(splus_prime) - sminus_prime * splus_prime * SQR(log(splus_prime))) / (128. * M_PI * sminus_prime * splus_prime);
          };
          alpha_pp *= std::norm(U[flav][k]);

          if (majorana) // For Majorana fermions, we can scatter off neutrinos and off antineutrinos
            alpha_pp *= 2;

          alpha_pp *= 2; // Each scattering generates at least 2 neutrinos
          if (majorana)  // For Majorana fermions, all final states are observable
            alpha_pp *= 2;
        }
        tot += SQR(SQR(mphi)) / (2 * mn[k]) * alpha_pp;

        if (alpha_s < 0 || alpha_t / SQR(SQR(g / mphi)) < -1e-11 || alpha_u / SQR(SQR(g / mphi)) < -1e-11 || alpha_tu / SQR(SQR(g / mphi)) < -1e-11 || (alpha_st + alpha_s + alpha_t) / SQR(SQR(g / mphi)) < -1e-11)
        {
          std::cerr << "Negative cross section when computing alpha; for tminus/mphi^2 = "
                    << tminus << ", tplus/mphi^2 = " << tplus << ", sminus_prime/mphi^2 = " << sminus_prime
                    << ", splus_prime/mphi^2 = " << splus_prime << ". The values of alpha are as follows:" << std::endl;
          std::cerr << "alpha_s = " << alpha_s << std::endl;
          std::cerr << "alpha_t = " << alpha_t << std::endl;
          std::cerr << "alpha_u = " << alpha_u << std::endl;
          std::cerr << "alpha_tu = " << alpha_tu << std::endl;
          std::cerr << "alpha_st + alpha_s + alpha_t = " << alpha_st + alpha_s + alpha_t << std::endl;
          std::cerr << "Possible roundoff errors for g=" << g << ", mphi=" << mphi << ", mntot=" << mntot << std::endl;
        }
      }

      return tot;
    }

  }; // End class calculate_flux

} // End namespace nuSIprop
#endif
