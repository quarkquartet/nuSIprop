#ifndef _INTERP_H_
#define _INTERP_H_

#include <iostream>
#include <cstring>
#include <cmath>
#include <unistd.h>
#define SQR(x)  ((x)*(x))  // square of a number
#define CUB(x)  ((x)*(x)*(x))  // cube of a number

namespace interp{

  template <int N_dim>
  class spline_ND{
    /**
     * This class implements N-dimensional spline interpolation. For reference, see
     *   https://en.wikipedia.org/wiki/Cubic_Hermite_spline
     *   https://en.wikipedia.org/wiki/Multivariate_interpolation
     *
     * Usage: The template parameter N_dim is the number of arguments of the function
     *        we want to interpolate.
     *        After passing the relevant variables to the constructor (see its
     *        documentation), just call f_eval(x0), where x0 is the N_dim-dimensional
     *        point at which we want to interpolate the function.
     *        For the moment, no extrapolation is implemented, and the code exits
     *        with an error if extrapolation is asked for.
     *
     * Example code:
     *
     *        // Declare arrays with nodes
     *        int N_dat[2] = {20, 18}; // 20 nodes in the x direction, 18 in the y direction
     *        double* data[2];
     *        data[0] = new double[N_dat[0]];
     *        data[1] = new double[N_dat[1]];
     *        double f[N_dat[0]][N_dat[1]];
     *
     *        // Fill nodes and f
     *        for(int i=0; i<20; ++i){
     *          data[0][i] = 2*M_PI * i/19.0;
     *          for(int j=0; j<18; ++j){
     *            data[1][j] = 2*M_PI * j/17.0;
     *            f[i][j] = cos(data[0][i] + data[1][j]);
     *          }
     *        }
     *
     *        // Create the interpolating object
     *        interp::spline_ND<2> spl(N_dat, data, f, true);
     *
     *        // Interpolate
     *        double x0[2] = {0.4, 1.6}
     *        interp = spl.f_eval(x0)
     *        
     * Author: Ivan Esteban
     */    
  public:
    /**
     * Constructor
     *
     * Mandatory parameters:
     *   int N_dat[N_dim] --- 1D array of size N_dim. N_dat[i] is the number of interpolating
     *                        nodes in the dimension i
     *   double *x[N_dim] --- 1D array of pointers. x[i] points to a *sorted* 1-D array with
     *                        the coordinates of the nodes in the dimension i
     *   double f[][]...  --- N_dim-dimensional array of doubles with the value of the
     *                        interpolating function at the nodes. More explicitly,
     *                           f[i][j][k]... = f(x[0][i], x[1][j], x[2][k], ...)
     *                        where f is the function we want to interpolate
     * Optional parameters:
     *   isRegular        --- Whether the nodes are equally spaced (true) or not (false).
     *                        If isLog==true, isRegular=true iff the nodes are equally spaced *in logarithmic space*
     *                        Default: false
     *   isLog[N_dim+1]   --- If isLog[i] is true, interpolate f(log[x_i]).
     *                        If isLog[N_dim] is true, interpolate log[f(...)]
     *                        Default: false
     */
    spline_ND(const int N_dat_[N_dim], const double *x_[N_dim], const void *f_, bool isRegular_ = false,
	      const bool *isLog_ = nullptr):
      isRegular(isRegular_){
      
      /* We first copy all arrays */
      
      N_dat = new int[N_dim];
      memcpy(N_dat, N_dat_, N_dim * sizeof(int));

      x = new double *[N_dim];
      for(int i=0; i<N_dim; ++i){
	x[i] = new double[N_dat[i]];
	memcpy(x[i], x_[i], N_dat[i] * sizeof(double));
      }

      isLog = new bool[N_dim+1];
      if (isLog_ == nullptr)
	for(int i=0; i<N_dim+1; ++i)
	  isLog[i] = false;
      else
	memcpy(isLog, isLog_, (N_dim+1) * sizeof(bool));
      
      // Internally, f will be stored as a 1-D array, so that for instance for N_dim=3
      //  f_[i][j][k] = f[k + j*N_dat[2] + i*N_dat[2]*Ndat[1]]
      // This mimicks how C stores multidimensional arrays, and simplifies memory copying
      // We do this to easily work with a variable-depth array
      int size_f = 1; // Number of elements in f
      for(int i=0; i<N_dim; ++i)
	size_f *= N_dat[i];
      f = new double[size_f];
      memcpy(f, f_, size_f * sizeof(double));

      /* Reparametrizations for logarithmic interpolation */

      // Reparametrize nodes 
      for(int i=0; i<N_dim; ++i)
	if(isLog[i])
	  for(int j=0; j<N_dat[i]; ++j){
	    if (x[i][j] < 0)
	      std::cerr<<"Error at interp: logarithmic interpolation was asked,"
		       <<" but the node #"<<j+1<<" of the coordinate #"<<i+1
		       <<" , "<<x[i][j]<<", is negative."<<std::endl;
	    x[i][j] = log(x[i][j]);
	  }
      
      // Reparametrize function values 
      if(isLog[N_dim])
	for(int j=0; j<size_f; ++j){
	  if (f[j] < 0)
	    std::cerr<<"Error at interp: logarithmic interpolation was asked,"
		     <<" but the function value at the node #"<<j
		     <<" , "<<f[j]<<", is negative."<<std::endl;
	  f[j] = log(f[j]);
	}
      

      /* Now, we compute all weights */
      computeWeights();
    }

    spline_ND(){ // Empty constructor for cython compatibility
      isRegular = true;
      isLog = new bool[N_dim+1](); // The () fills isLog with false      
      N_dat = new int[N_dim](); // The () fills N_dat with zeros
      x = new double *[N_dim];
      for(int i=0; i<N_dim; ++i)
	x[i] = new double[0];

      f = new double[0];

      computeWeights();
    }
      
    /**
     * Alternative constructor. Reads the data directly from a file.
     *
     * Mandatory parameters:
     *   int N_dat[N_dim] --- 1D array of size N_dim. N_dat[i] is the number of interpolating
     *                        nodes in the dimension i
     *   char* filepath   --- Path of the file with the data. The syntax must be
     *                        x[0][0]        x[1][0]        ... x[N_dim-1][0] f
     *                        x[0][0]        x[1][0]        ... x[N_dim-1][1] f
     *                        x[0][0]        x[1][0]        ... x[N_dim-1][2] f
     *                        ...
     *                        x[0][N_dat[0]] x[1][N_dat[1]] ... x[N_dim-1][N_dat[N_dim-1]] f
     *
     * Optional parameters:
     *   isRegular        --- Whether the nodes are equally spaced (true) or not (false).
     *                        If isLog==true, isRegular=true iff the nodes are equally spaced *in logarithmic space*
     *                        Default: false
     *   isLog[N_dim+1]   --- If isLog[i] is true, interpolate f(log[x_i]).
     *                        If isLog[N_dim] is true, interpolate log[f(...)]
     *                        Default: false
     *   isBinary         --- Whether data is stored in binary format. In such case, the order is assumed to be
     *                        as for the text file, and the data is assumed to be *float*
     *                        Default: false
     */
    spline_ND(const int N_dat_[N_dim], const char *filepath, bool isRegular_ = false,
	      const bool *isLog_ = nullptr, bool isBinary = false):
      isRegular(isRegular_){
      
      /* We first copy the array N_dat and initialize other arrays */
      
      N_dat = new int[N_dim];
      memcpy(N_dat, N_dat_, N_dim * sizeof(int));

      x = new double *[N_dim];
      for(int i=0; i<N_dim; ++i)
	x[i] = new double[N_dat[i]];

      isLog = new bool[N_dim+1];
      if (isLog_ == nullptr)
	for(int i=0; i<N_dim+1; ++i)
	  isLog[i] = false;
      else
	memcpy(isLog, isLog_, (N_dim+1) * sizeof(bool));

      int size_f = 1; // Number of elements in f
      for(int i=0; i<N_dim; ++i)
	size_f *= N_dat[i];
      f = new double[size_f];

      /* We now read the input file */
      if(!isBinary){
	/* Text file */
	if(access(filepath, F_OK) != 0){
	  std::cerr<<"Error at interp: the input file "<<filepath<<" does not exist"<<std::endl;
	  exit(1);
	}
	
	FILE *fp = fopen(filepath, "r");
	int idx[N_dim] = {}; // This array will keep track of the index of each variable we are looping over.
	int p = N_dim-1; // p will keep track of which variable we are looping over
	char line[1000]; // Line to be scanned
	char *var; // Each variable to be read
	while(idx[0] < N_dat[0]){
	  fgets(line, sizeof line, fp);
	  if (line[0] == '#') continue; // Ignore lines starting with #
	  var = strtok(line, " "); // Read first element in line
	
	  for(int i=0; i<N_dim; ++i){
	    sscanf(var, "%lf", &x[i][idx[i]]);
	    var = strtok(NULL, " "); // Read next element in line
	  }

	  // Scan f
	  int idx_f = 0; // idx_f = idx[N_dim-1]
	  //       + idx[N_dim-2]*N_dat[N_dim-1]
	  //       + idx[N_dim-3]*N_dat[N_dim-1]*N_dat[N_dim-2]
	  //       + ...
	  for(int i=N_dim-1; i>=0; --i){
	    int tmp = idx[i];
	    for (int j=i+1; j<=N_dim-1; ++j)
	      tmp *= N_dat[j];
	    idx_f += tmp;
	  }
	  sscanf(var, "%lf", &f[idx_f]);
	  var = strtok(NULL, " "); // Read next element in line

	  // Advance the loop               
          ++idx[p];
          while(idx[p]==N_dat[p]){ // Iteration over p is over. Look for the next index that has to be iterated over
            if(p==0)
              break;
            idx[p] = 0; // Iteration over p is over, reset it to 0
            p -= 1; // Decrease p by 1         
            ++idx[p]; // Increase the p+1-th index
          }
          if(idx[0] == N_dat[0])
            break;
          p = N_dim-1;
	}
	fclose(fp);
      } else{
	/* Binary file */
	if(access(filepath, F_OK) != 0){
	  std::cerr<<"Error at interp: the input file "<<filepath<<" does not exist"<<std::endl;
	  exit(1);
	}
	FILE *fp = fopen(filepath, "rb");
	int idx[N_dim] = {}; // This array will keep track of the index of each variable we are looping over.
	int p = N_dim-1; // p will keep track of which variable we are looping over
	float var[N_dim+1]; //Variables to be read
	while(idx[0] < N_dat[0]){
	  fread(&var, sizeof(float), N_dim+1, fp);

	  for(int i=0; i<N_dim; ++i)
	    x[i][idx[i]] = (double) var[i];
	  
	  // Scan f
	  int idx_f = 0; // idx_f = idx[N_dim-1]
	  //       + idx[N_dim-2]*N_dat[N_dim-1]
	  //       + idx[N_dim-3]*N_dat[N_dim-1]*N_dat[N_dim-2]
	  //       + ...
	  for(int i=N_dim-1; i>=0; --i){
	    int tmp = idx[i];
	    for (int j=i+1; j<=N_dim-1; ++j)
	      tmp *= N_dat[j];
	    idx_f += tmp;
	  }
	  f[idx_f] = (double) var[N_dim];

	  // Advance the loop               
          ++idx[p];
          while(idx[p]==N_dat[p]){ // Iteration over p is over. Look for the next index that has to be iterated over
            if(p==0)
              break;
            idx[p] = 0; // Iteration over p is over, reset it to 0
            p -= 1; // Decrease p by 1         
            ++idx[p]; // Increase the p+1-th index
          }
          if(idx[0] == N_dat[0])
            break;
          p = N_dim-1;
	}
	fclose(fp);
      }
      
      /* Reparametrizations for logarithmic interpolation */

      // Reparametrize nodes
      for(int i=0; i<N_dim; ++i)
	if(isLog[i])
	  for(int j=0; j<N_dat[i]; ++j){
	    if (x[i][j] < 0)
	      std::cerr<<"Error at interp: logarithmic interpolation was asked,"
		       <<" but the node #"<<j+1<<" of the coordinate #"<<i+1
		       <<" , "<<x[i][j]<<", is negative."<<std::endl;
	    x[i][j] = log(x[i][j]);
	  }
      
      // Reparametrize function values 
      if(isLog[N_dim])
	for(int j=0; j<size_f; ++j){
	  if (f[j] < 0)
	    std::cerr<<"Error at interp: logarithmic interpolation was asked,"
		     <<" but the function value at the node #"<<j
		     <<" , "<<f[j]<<", is negative."<<std::endl;
	  f[j] = log(f[j]);
	}
      

      /* Now, we compute all weights */
      computeWeights();
    }

    /**
     * Brace initializator versions of constructors
     */

    spline_ND(const std::initializer_list<int>& N_dat_, const double *x_[N_dim], const void *f_,
	      bool isRegular_ = false, const std::initializer_list<bool>& isLog_ = {}) :
      spline_ND(N_dat_.begin(), x_, f_, isRegular_, (isLog_.size()==0)? nullptr : isLog_.begin()) {}
    
    spline_ND(const std::initializer_list<int>& N_dat_, const char *filepath, bool isRegular_ = false,
	      const std::initializer_list<bool>& isLog_ = {},
	      bool isBinary = false) :
      spline_ND(N_dat_.begin(), filepath, isRegular_, (isLog_.size()==0)? nullptr : isLog_.begin(),
		isBinary) {}

    /**
     * Returns the interpolated value for any point x0
     *
     *   double x0[N_dim] --- Point at which we want to interpolate the function
     *
     *   For the moment, no extrapolation is implemented, and the code exits
     *   with an error if extrapolation is asked for.
     *   I.e., we must have x[i][0] < x0[i] < x[i][N_dat[i]-1]
     */
    double f_eval(const double x0_[N_dim]) const{
      double *x0 = new double[N_dim];
      memcpy(x0, x0_, N_dim * sizeof(double));

      /* Reparametrize if there is logarithmic interpolation */
      for(int i=0; i<N_dim; ++i)
	if(isLog[i])
	  x0[i] = log(x0[i]);
      
      /* Check that x0 is within bounds */
      for(int i=0; i<N_dim; ++i) 
	if ((x0[i] <= x[i][0]) || (x0[i] >= x[i][N_dat[i]-1])){
	    std::cerr<<"Error at interp: the coordinate #"<<i+1
		     <<" of the input value, "<<x0[i]<<", is outside its bounds ("
		     <<x[i][0]<<", "<<x[i][N_dat[i]-1]<<")"<<std::endl;
	    exit(1);
	  }

      /* Find the position of x0 in the interpolating nodes */
      int k[N_dim]; // x[i][k[i]] < x0[i] < x[i][k[i]+1]

      if(isRegular) // Nodes are equally spaced: x[i][k] = x[i][0] + (x[i][1]-x[i][0]) * k
	for(int i=0; i<N_dim; ++i){
	  k[i] = (x0[i] - x[i][0]) / (x[i][1]-x[i][0]);
	  // Set extreme cases by hand to avoid roundoff uncertainties
	  if(x0[i] < x[i][1])
	    k[i] = 0;
	  else if(x0[i] > x[i][N_dat[i]-2])
	    k[i] = N_dat[i]-2;
	}
      else 
	for(int i=0; i<N_dim; ++i){
	  // Perform a binary search
	  int L = 0, R = N_dat[i] - 1;
	  int m;
	  while (L <= R){
	    m =  (L+R) / 2;
	    if (x0[i] < x[i][m])
	      R = m - 1;
	    else{
	      k[i] = m;
	      L = m + 1;
	    }
	  }
	}

      /* Obtain the nodes at which the function will be evaluated */
      int idx_min[N_dim+1], idx_max[N_dim+1]; 
      
      for(int i=0; i<N_dim; ++i){
	if (k[i] == 0) {
	  idx_min[i] = k[i];
	  idx_max[i] = k[i]+2;
	} else if (k[i] == N_dat[i]-2) {
	  idx_min[i] = k[i]-1;
	  idx_max[i] = k[i]+1;
	} else{
	  idx_min[i] = k[i]-1;
	  idx_max[i] = k[i]+2;
	}
      }

      // The extra element is necessary for the N-D loop that computes the output
      idx_min[N_dim] = 0;
      idx_max[N_dim] = 1;
      
      /* Interpolate */
      double t[N_dim];
      for(int i=0; i<N_dim; ++i)
	t[i] = (x0[i] - x[i][k[i]]) / (x[i][k[i]+1] - x[i][k[i]]);

      /*
       * We want to compute
       * ( \sum_{idx_0=0}^{idx_max[0]-idx_min[0]} ) ( \sum_{idx_1=0}^{idx_max[1]-idx_min[1]} ) (...)
       *    f[idx_min[0] + idx_0][idx_min[1] + idx_1][...] *
       *       (t[0]^3 * w[0][idx][0][k[0]] + t[0]^2 * w[0][idx][1][k[0]] + t[0] * w[0][idx][2][k[0]] + w[0][idx][3][k[0]]) * 
       *       (t[1]^3 * w[1][idx][0][k[1]] + t[1]^2 * w[1][idx][1][k[1]] + t[1] * w[1][idx][2][k[1]] + w[1][idx][3][k[1]]) *
       *       (...)
       */

      double res = 0; // Result
      int idx[N_dim+1]; // This array will keep track of the index of each variable we are looping over.
			// The extra element will be 1 when all iterations are over.
      for(int i=0; i<N_dim+1; ++i)
	idx[i] = 0;
      
      int p = 0; // p will keep track of which variable we are looping over
      while(idx[N_dim]==0){
	// Do the work 
	int idx_f = 0; // idx_f = (idx_min[N_dim-1] + idx[N_dim-1])
		       //       + (idx_min[N_dim-2] + idx[N_dim-2])*N_dat[N_dim-1]
		       //       + (idx_min[N_dim-3] + idx[N_dim-3])*N_dat[N_dim-1]*N_dat[N_dim-2]
		       //       + ...
	for(int i=N_dim-1; i>=0; --i){
	  int tmp = idx_min[i] + idx[i];
	  for (int j=i+1; j<=N_dim-1; ++j)
	    tmp *= N_dat[j];
	  idx_f += tmp;
	}
	
	double tmp = f[idx_f];
	for(int i=0; i<N_dim; ++i)
	  tmp *= (CUB(t[i]) * w[i][idx[i]][0][k[i]] + SQR(t[i]) * w[i][idx[i]][1][k[i]]
		  + t[i] * w[i][idx[i]][2][k[i]] + w[i][idx[i]][3][k[i]]);
	res += tmp;

	// Advance the loop
	++idx[0];
	while(idx[p]==idx_max[p]-idx_min[p]+1){ // Iteration over p is over. Look for the next index that has to be iterated over
	  idx[p] = 0; // Iteration over p is over, reset it to 0
	  p += 1; // Increase p by 1
	  ++idx[p]; // Increase the p+1-th index
	  if(idx[p]!=idx_max[p]-idx_min[p]+1) // If we are not done with iterating the p+1-th index, we have to start over to iterate it
	    p = 0;
	}
      }

      delete[] x0;
      if(isLog[N_dim])
	return exp(res);
      else
	return res;
    }

    /**
     * Brace initializator version of f_eval
     */

    double f_eval(const std::initializer_list<double>& x0_) const{
      if (x0_.size() != N_dim){
	std::cerr<<"Error at interp: a function of dimension"<<N_dim
		 <<" was asked to be interpolated at a point with dimension"
		 <<x0_.size()<<std::endl;
	exit(1);
      }
      return f_eval(x0_.begin());
    }

    // Copy constructor
    spline_ND(const spline_ND<N_dim> &orig):
      isRegular(orig.isRegular){
      /* Properly take care of memory in the arrays */
      N_dat = new int[N_dim];
      memcpy(N_dat, orig.N_dat, N_dim * sizeof(int));

      x = new double *[N_dim];
      for(int i=0; i<N_dim; ++i){
	x[i] = new double[N_dat[i]];
	memcpy(x[i], orig.x[i], N_dat[i] * sizeof(double));
      }

      isLog = new bool[N_dim+1];
      memcpy(isLog, orig.isLog, (N_dim+1) * sizeof(bool));

      int size_f = 1; // Number of elements in f
      for(int i=0; i<N_dim; ++i)
	size_f *= N_dat[i];
      f = new double[size_f];
      memcpy(f, orig.f, size_f * sizeof(double));

      /* Compute all weights */
      computeWeights();
    }

    // Assignment operator
    spline_ND<N_dim>& operator=(const spline_ND<N_dim> &rhs){
      isRegular = rhs.isRegular;
      
      /* Properly take care of memory in the arrays */
      delete[] N_dat;
      for(int i=0; i<N_dim; ++i)
	delete[] x[i];
      delete[] x;
      delete[] isLog;
      delete[] f;
      for(int i=0; i<N_dim; ++i)
	for(int j=0; j<4; ++j)
	  for(int k=0; k<4; ++k)
	    delete[] w[i][j][k];
      
      N_dat = new int[N_dim];
      memcpy(N_dat, rhs.N_dat, N_dim * sizeof(int));

      x = new double *[N_dim];
      for(int i=0; i<N_dim; ++i){
	x[i] = new double[N_dat[i]];
	memcpy(x[i], rhs.x[i], N_dat[i] * sizeof(double));
      }

      isLog = new bool[N_dim+1];
      memcpy(isLog, rhs.isLog, (N_dim+1) * sizeof(bool));

      int size_f = 1; // Number of elements in f
      for(int i=0; i<N_dim; ++i)
	size_f *= N_dat[i];
      f = new double[size_f];
      memcpy(f, rhs.f, size_f * sizeof(double));

      /* Compute all weights */
      computeWeights();
      return *this;
    }

    // Destructor
    ~spline_ND(){
      delete[] N_dat;
      
      for(int i=0; i<N_dim; ++i)
	delete[] x[i];
      delete[] x;

      delete[] isLog;

      delete[] f;

      for(int i=0; i<N_dim; ++i)
	for(int j=0; j<4; ++j)
	  for(int k=0; k<4; ++k)
	    delete[] w[i][j][k];
    }
  private:
    int *N_dat; // Number of interpolating nodes in the dimension i
    double **x; // x[i][j] is the coordinate of the j-th node in the dimension i
    double *f; // Function values at nodes, stored as a flattened 1D array (see the constructor)
    double *w[N_dim][4][4]; // Interpolating weights
    bool isRegular; // Equally spaced nodes?
    bool *isLog; // Logarithmic interpolation?

    /**
     * This function computes the interpolating weights
     */
    void computeWeights(){
      for(int i=0; i<N_dim; ++i)
	for(int j=0; j<4; ++j)
	  for(int k=0; k<4; ++k)
	    w[i][j][k] = new double[N_dat[i]];

      for(int i=0; i<N_dim; ++i)
	for(int j=0; j<N_dat[i]; ++j){
	  if (j==0){
	    w[i][0][0][j] = 0;
	    w[i][0][1][j] = (x[i][j]-x[i][j+1]) / (x[i][j]-x[i][j+2]);
	    w[i][0][2][j] = (-1 + (x[i][j+1]-x[i][j]) / (x[i][j]-x[i][j+2]));
	    w[i][0][3][j] = 1;

	    w[i][1][0][j] = 0;
	    w[i][1][1][j] = (x[i][j+1]-x[i][j]) / (x[i][j+1]-x[i][j+2]);
	    w[i][1][2][j] = (x[i][j]-x[i][j+2]) / (x[i][j+1]-x[i][j+2]);
	    w[i][1][3][j] = 0;

	    w[i][2][0][j] = 0;
	    w[i][2][1][j] = SQR(x[i][j+1]-x[i][j]) / ((x[i][j+2]-x[i][j+1])*(x[i][j+2]-x[i][j]));
	    w[i][2][2][j] = SQR(x[i][j+1]-x[i][j]) / ((x[i][j+2]-x[i][j+1])*(x[i][j]-x[i][j+2]));
	    w[i][2][3][j] = 0;
	  } else if (j==N_dat[i]-2){
	    w[i][0][0][j] = 0;
	    w[i][0][1][j] = SQR(x[i][j+1]-x[i][j]) / ((x[i][j-1]-x[i][j])*(x[i][j-1]-x[i][j+1]));
	    w[i][0][2][j] = SQR(x[i][j+1]-x[i][j]) / ((x[i][j]-x[i][j-1])*(x[i][j-1]-x[i][j+1]));
	    w[i][0][3][j] = 0;

	    w[i][1][0][j] = 0;
	    w[i][1][1][j] = (x[i][j+1]-x[i][j])/(x[i][j-1]-x[i][j]);
	    w[i][1][2][j] = (2*x[i][j]-x[i][j+1]-x[i][j-1])/(x[i][j-1]-x[i][j]);
	    w[i][1][3][j] = 1;

	    w[i][2][0][j] = 0;
	    w[i][2][1][j] = (x[i][j]-x[i][j+1])/(x[i][j-1]-x[i][j+1]);
	    w[i][2][2][j] = (x[i][j-1]-x[i][j])/(x[i][j-1]-x[i][j+1]);
	    w[i][2][3][j] = 0;
	  } else{
	    w[i][0][0][j] = SQR(x[i][j+1]-x[i][j])/((x[i][j] - x[i][j-1])*(x[i][j-1] - x[i][j+1]));
	    w[i][0][1][j] = 2*SQR(x[i][j+1]-x[i][j])/((x[i][j-1] - x[i][j])*(x[i][j-1] - x[i][j+1]));
	    w[i][0][2][j] = SQR(x[i][j+1]-x[i][j])/((x[i][j] - x[i][j-1])*(x[i][j-1] - x[i][j+1]));
	    w[i][0][3][j] = 0;
	    
	    w[i][1][0][j] = (x[i][j]-x[i][j+1]) * (1/(x[i][j-1]-x[i][j]) + 1/(x[i][j]-x[i][j+2]));
	    w[i][1][1][j] = (x[i][j]-x[i][j+1]) * (2/(x[i][j]-x[i][j-1]) + 1/(x[i][j+2]-x[i][j]));
	    w[i][1][2][j] = (2*x[i][j]-x[i][j+1]-x[i][j-1])/(x[i][j-1]-x[i][j]);
	    w[i][1][3][j] = 1;

	    w[i][2][0][j] = (x[i][j+1]-x[i][j]) * (1/(x[i][j-1]-x[i][j+1]) + 1/(x[i][j+1]-x[i][j+2]));
	    w[i][2][1][j] = (x[i][j+1]-x[i][j]) * (2/(x[i][j+1]-x[i][j-1]) + 1/(x[i][j+2]-x[i][j+1]));
	    w[i][2][2][j] = (x[i][j-1]-x[i][j])/(x[i][j-1]-x[i][j+1]);
	    w[i][2][3][j] = 0;

	    w[i][3][0][j] = SQR(x[i][j+1]-x[i][j])/((-x[i][j+1] + x[i][j+2])*(-x[i][j] + x[i][j+2]));
	    w[i][3][1][j] = SQR(x[i][j+1]-x[i][j])/((x[i][j+1] - x[i][j+2])*(-x[i][j] + x[i][j+2]));
	    w[i][3][2][j] = 0;
	    w[i][3][3][j] = 0;
	  }	  
	}
    }
    
  }; // End class spline_ND

} // End namespace interp

#endif
