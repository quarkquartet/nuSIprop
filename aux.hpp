#include <gsl/gsl_poly.h>
#include <gsl/gsl_sf_dilog.h>
#include <complex.h>

/**
 * Auxiliary functions and constants for nuSIprop
 * Author: Ivan Esteban
 */

namespace nuSIaux{

  double getmL(double mSum, double dmqSL, double dmqAT){
    /** Returns the smallest neutrino mass
     *  mSum  --- Total neutrino mass
     *  dmqSL --- \delta m^2_{21}
     *  dmqAT --- \delta m^2_{31} for NO, \delta m^2_{32} for IO
     */
    double coeff[5];
    if(dmqAT>0){
      coeff[0] = -( SQR(SQR(mSum))-2*SQR(mSum)*(dmqSL+dmqAT)+SQR(dmqSL-dmqAT) );
      coeff[1] = 4*mSum*(SQR(mSum)-dmqSL-dmqAT);
      coeff[2] = -2*(SQR(mSum)-dmqSL-dmqAT);
      coeff[3] = -4*mSum;
      coeff[4] = 3;
    } else{
      coeff[0] = -( SQR(SQR(mSum))-2*SQR(mSum)*(-dmqSL-2*dmqAT)+SQR(dmqSL) );
      coeff[1] = 4*mSum*(SQR(mSum)+dmqSL+2*dmqAT);
      coeff[2] = -2*(SQR(mSum)+dmqSL+2*dmqAT);
      coeff[3] = -4*mSum;
      coeff[4] = 3;
    }
    double sol[10];

    gsl_poly_complex_workspace * w = gsl_poly_complex_workspace_alloc(5);
    gsl_poly_complex_solve(coeff, 5, w, sol);
    gsl_poly_complex_workspace_free(w);
  
    for(int i=0;i<4;++i){
      double ml = sol[2*i];
      if(	(sol[2*i+1]<1.0e-7) && (ml >= 0) && //Real positive solution
		(mSum-ml>1.0e-7) && //1st constraint
		( (dmqAT>0 &&(SQR(mSum)-dmqAT-dmqSL-SQR(ml)-2*ml*mSum > 1.0e-7) ) ||
		  (dmqAT<0 &&(SQR(mSum)+2*dmqAT+dmqSL-SQR(ml)-2*ml*mSum > 1.0e-7) ) )//2nd constraint
		)
	return ml;
    }

    std::cerr<<"No neutrino mass spectrum was found corresponding to \\sum m = "<<mSum<<", dmqAT = "<<dmqAT<<", dmqSL = "<<dmqSL<<". Exiting..."<<std::endl;
    exit(1);
  }

  /* Weights and nodes for 3-point Gauss-Legendre quadrature */
  double w_integ[3] = {5./9., 8./9., 5./9.};
  double x_integ[3] = {-sqrt(3./5.), 0, sqrt(3./5.)};

  /* Fast and precise implementations of f(x)-f(y) for various f
     The logic is that there are 2 main improvements wrt a naive computation of the integrals in nuSIprop:
      1-Improve the speed by approximating dilogarithms/atan by Taylor expansions when arguments are large/small
      2-Improve numerical convergence by using Taylor expansions for (1+small argument)
     these are implemented below
   */

  double atandiff(double x, double y){
    /**
     * Fast implementation of atan(x)-atan(y) if x, y >> 1.
     * The relative error we introduce here by truncating the Taylor
     * series of atan is O(10^-10), probably larger than the precision
     * error if we computed atan(x)-atan(y)
     */
    if (fabs(x)<1e2 || fabs(y)<1e2 || x*y<0)
      return atan(x) - atan(y);
    else
      return -1./x + CUB(1./x)/3. -
	(-1./y + CUB(1./y)/3.);
  }

  double _Complex dilogdiff_complex(double _Complex x, double _Complex y){
    /**
     * Fast, precise implementation of dilog(x)-dilog(y) if |x|, |y| >> 1 for complex x, y
     * The relative error we introduce here by truncating the Taylor
     * series of dilog is O(10^-10), probably larger than the precision
     * error if we computed dilog(x)-dilog(y)
     */
    if (cabs(x) > 1e2 && cabs(y) > 1e2){
      int sign_im_x = (cimag(x) >= 0)? 1 : -1;
      int sign_im_y = (cimag(y) >= 0)? 1 : -1;      
      return -1/(16.*SQR(SQR(x))) - 1/(9.*CUB(x)) - 1/(4.*SQR(x)) - 1/x - I/2*(-sign_im_x*2*M_PI*clog(x) - I*SQR(clog(x))) -
	(-1/(16.*SQR(SQR(y))) - 1/(9.*CUB(y)) - 1/(4.*SQR(y)) - 1/y - I/2*(-sign_im_y*2*M_PI*clog(y) - I*SQR(clog(y))));
    }
    else{
      gsl_sf_result dilog_x_re, dilog_x_im, dilog_y_re, dilog_y_im;
      gsl_sf_complex_dilog_xy_e(creal(x), cimag(x), &dilog_x_re, &dilog_x_im);
      gsl_sf_complex_dilog_xy_e(creal(y), cimag(y), &dilog_y_re, &dilog_y_im);
      return dilog_x_re.val + I*dilog_x_im.val - dilog_y_re.val - I*dilog_y_im.val;
    }
  }

  double dilogdiff(double x, double y){
    /**
     * Fast, precise implementation of dilog(-x)-dilog(-y) if x, y >> 1 or x, y << 1.
     * The relative error we introduce here by truncating the Taylor
     * series of dilog is O(10^-10), probably larger than the precision
     * error if we computed dilog(x)-dilog(y)
     */
    if (x > 1e2 && y > 1e2)
      return -SQR(log(x))/2. + 1./x - SQR(1./x)/4. + CUB(1./x)/9. - SQR(SQR(1./x))/16 -
	(-SQR(log(y))/2. + 1./y - SQR(1./y)/4. + CUB(1./y)/9. - SQR(SQR(1./y))/16);
    else if (x < 1e-2 && y < 1e-2)
      return -x + SQR(x)/4. - CUB(x)/9. + SQR(SQR(x))/16. -
	(-y + SQR(y)/4. - CUB(y)/9. + SQR(SQR(y))/16.);
    else
      return gsl_sf_dilog(-x) - gsl_sf_dilog(-y);
  }

  double dilog1mdiff(double x, double y){
    /**
     * Fast, precise implementation of dilog(-1-x)-dilog(-1-y) if x, y >> 1 or x, y << 1.
     * The relative error we introduce here by truncating the Taylor
     * series of dilog is O(10^-10), probably larger than the precision
     * error if we computed dilog(-1-x)-dilog(-1-y)
     */
    if (x > 1e2 && y > 1e2)
      return -SQR(log(x))/2. + (1 - log(x))/x + (-7 + 2*log(x))/(4.*SQR(x)) + (19 - 3*log(x))/(9.*CUB(x)) + (-125 + 12*log(x))/(48.*SQR(SQR(x))) -
	(-SQR(log(y))/2. + (1 - log(y))/y + (-7 + 2*log(y))/(4.*SQR(y)) + (19 - 3*log(y))/(9.*CUB(y)) + (-125 + 12*log(y))/(48.*SQR(SQR(y))));
    else if (x < 1e-2 && y < 1e-2)
      return -x*log(2) + (SQR(x)*(-1 + log(4)))/4. + (CUB(x)*(5 - 8*log(2)))/24. + SQR(SQR(x))*(-1./6. + log(2)/4.) -
	(-y*log(2) + (SQR(y)*(-1 + log(4)))/4. + (CUB(y)*(5 - 8*log(2)))/24. + SQR(SQR(y))*(-1./6. + log(2)/4.));
    else
      return gsl_sf_dilog(-1-x) - gsl_sf_dilog(-1-y);
  }

  double dilog1pdiff(double x, double y){
    /**
     * Fast, precise implementation of dilog(1+x)-dilog(1+y) if |x|, |y| >> 1 or |x|, |y| << 1.
     * It is assumed x<0, y<0.
     * The relative error we introduce here by truncating the Taylor
     * series of dilog is O(10^-10), probably larger than the precision
     * error if we computed dilog(1+x)-dilog(1+y)
     */
    if (-x > 1e2 && -y > 1e2)
      return (-1 - 3*log(-x))/(9.*CUB(x)) + (-1 - log(-x))/x - SQR(log(-x))/2. + (1 + 2*log(-x))/(4.*SQR(x)) + (1 + 4*log(-x))/(16.*SQR(SQR(x))) -
	((-1 - 3*log(-y))/(9.*CUB(y)) + (-1 - log(-y))/y - SQR(log(-y))/2. + (1 + 2*log(-y))/(4.*SQR(y)) + (1 + 4*log(-y))/(16.*SQR(SQR(y))));
    else if (-x < 1e-2 && -y < 1e-2)
      return x*(1 - log(-x)) + (SQR(x)*(-1 + 2*log(-x)))/4. + (CUB(x)*(1 - 3*log(-x)))/9. + (SQR(SQR(x))*(-1 + 4*log(-x)))/16. -
	(y*(1 - log(-y)) + (SQR(y)*(-1 + 2*log(-y)))/4. + (CUB(y)*(1 - 3*log(-y)))/9. + (SQR(SQR(y))*(-1 + 4*log(-y)))/16.);
    else
      return gsl_sf_dilog(1+x) - gsl_sf_dilog(1+y);
  }

  double dilog1over1mdiff(double x, double y){
    /**
     * Fast, precise implementation of dilog(1/(1-x))-dilog(1/(1-y)) if |x|, |y| >> 1 or |x|, |y| << 1.
     * It is assumed x<0, y<0.
     * The relative error we introduce here by truncating the Taylor
     * series of dilog is O(10^-10), probably larger than the precision
     * error if we computed dilog(1/(1-x))-dilog(1/(1-y))
     */
    if (-x > 1e2 && -y > 1e2)
      return -25/(48.*SQR(SQR(x))) - 11/(18.*CUB(x)) - 3/(4.*SQR(x)) - 1/x -
	(-25/(48.*SQR(SQR(y))) - 11/(18.*CUB(y)) - 3/(4.*SQR(y)) - 1/y);
    else if (-x < 1e-2 && -y < 1e-2)
      return (SQR(SQR(x))*(-19 - 12*log(-x)))/48. + (CUB(x)*(-7 - 6*log(-x)))/18. + (SQR(x)*(-1 - 2*log(-x)))/4. + x*(1 - log(-x)) -
	((SQR(SQR(y))*(-19 - 12*log(-y)))/48. + (CUB(y)*(-7 - 6*log(-y)))/18. + (SQR(y)*(-1 - 2*log(-y)))/4. + y*(1 - log(-y)));
    else
      return gsl_sf_dilog(1/(1-x)) - gsl_sf_dilog(1/(1-y));
  }

} // End namespace nuSIaux
