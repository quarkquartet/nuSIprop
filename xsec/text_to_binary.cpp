#include <fstream>
#include <cstring>

// Convert text files to binary files

int main(){
  FILE *fp = fopen("alpha_phiphi.dat", "r");

  char line[1000]; // Line to be scanned
  char *var; // Each variable to be read
  constexpr int n_lines = 100000000; // Number of lines to be read [excluding comments]
  constexpr int n_els = 4;
  float** data = new float*[n_lines];
  for(int i=0; i<n_lines; ++i)
    data[i] = new float[n_els];

  /* We read the data file */
  for(int j=0; j<n_lines; ++j){
    fgets(line, sizeof line, fp);
    if (line[0] == '#'){ // Ignore lines starting with #
      j -= 1;
      continue; 
    }
    
    var = strtok(line, " "); // Read first element in line
	
    for(int i=0; i<n_els; ++i){
      sscanf(var, "%f", &data[j][i]);
      var = strtok(NULL, " "); // Read next element in line
    }
  }
  fclose(fp);

  /* We write it to a binary file */
  FILE *f = fopen("alpha_phiphi.bin", "wb");
  for(int j=0; j<n_lines; ++j)
    fwrite(data[j], sizeof(float), n_els, f);
  fclose(f);

  for(int i=0; i<n_lines; ++i)
    delete[] data[i];
  delete[] data;

  /* And repeat */

  fp = fopen("alphatilde_phiphi.dat", "r");

  constexpr int n_lines_tilde = 500000; // Number of lines to be read [excluding comments]
  constexpr int n_els_tilde = 3;
  data = new float*[n_lines_tilde];
  for(int i=0; i<n_lines_tilde; ++i)
    data[i] = new float[n_els_tilde];

  /* We read the data file */
  for(int j=0; j<n_lines_tilde; ++j){
    fgets(line, sizeof line, fp);
    if (line[0] == '#'){ // Ignore lines starting with #
      j -= 1;
      continue; 
    }
    
    var = strtok(line, " "); // Read first element in line
	
    for(int i=0; i<n_els_tilde; ++i){
      sscanf(var, "%f", &data[j][i]);
      var = strtok(NULL, " "); // Read next element in line
    }
  }
  fclose(fp);

  /* We write it to a binary file */
  f = fopen("alphatilde_phiphi.bin", "wb");
  for(int j=0; j<n_lines_tilde; ++j)
    fwrite(data[j], sizeof(float), n_els_tilde, f);
  fclose(f);

  for(int i=0; i<n_lines_tilde; ++i)
    delete[] data[i];
  delete[] data;
}
