/* In this file, we define in C some functions to speed up computation
  Compile it as
  
  gcc -shared -o funcs.so -fPIC funcs.c
  */

#include <math.h>
#include <stdio.h>
#define SQR(x) ((x)*(x))
#define CUB(x) ((x)*(x)*(x))

double primitive(double tau, double s){
  /**
   * \int dtau/(-tau) dsigma/dtau
   */
  return (1/(1 + tau) + 1/((-1 + s)*(-1 + s + tau)) +
	  (-(SQR(-1 + s)*(4 + (-3 + s)*s)*log(-1 - tau)) + (-2 + s)*CUB(s)*log(-tau) + (-4 + s*(9 + (-5 + s)*s))*log(-1 + s + tau))/((-2 + s)*SQR(-1 + s)))
    /(64.*M_PI*SQR(s));
}

double dsigma_phiphi_over_tauphi(int n, double args[n]){
  double tbar = args[1];
  double sbar = args[0];

  /* if (sbar <= 4) */
  /*   return 0; */
  /* if (tauphibar <= -1 - 1./4.*SQR(sqrt(sbar) + sqrt(sbar-4))) */
  /*   return 0; */
  /* if (tauphibar >= -1 - 1./4.*SQR(sqrt(sbar) - sqrt(sbar-4))) */
  /*   return 0; */
  double upperLim = (tbar < -1 - 1./4.*SQR(sqrt(sbar) - sqrt(sbar-4)))? tbar : -1 - 1./4.*SQR(sqrt(sbar) - sqrt(sbar-4));
  double lowerLim = -1 - 1./4.*SQR(sqrt(sbar) + sqrt(sbar-4));
  if (upperLim < lowerLim){
    //printf("%.3e %.3e %.3e %.3e %.3e %.3e\n", tbar, sbar, upperLim, lowerLim, -SQR(tbar)/(1+tbar), upperLim-lowerLim);
    return 0;
  }
  
  return primitive(upperLim, sbar) - primitive(lowerLim, sbar);
}
