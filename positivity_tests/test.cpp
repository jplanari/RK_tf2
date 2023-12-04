#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>

#include "../mods/fsm_rk/butcherTableaux.h"

double f(double phi)
{
  return -1.0*phi;
}

double advanceRungeKutta(double phi, butcherTableau method, double h)
{
    double phin;
    int s = method.b.size();
    std::vector<double> phis(s);
    
    for(int i=0; i<s; ++i){
      phis.at(i) = phi;
      for(int j=0;j<i;++j)
        phis.at(i) += h*method.A.at(id(i+1,j+1))*f(phis.at(j));
    }
    
    phin = phi;
    for(int i=0; i<s; ++i)
      phin += h*method.b.at(i)*f(phis.at(i));

  return phin;
}


int main(int argc, char** argv){
  std::string schemeName = argv[1];
  butcherTableau method = intScheme(schemeName);
  displayMethod(method);
  double h = atof(argv[2]);
  double phi0 = 1.0;
  double phi=phi0;
  int N = atoi(argv[3]);
  std::string filename = schemeName+"/testPositivity_"+schemeName+"_h_"+argv[2]+"_N_"+argv[3]+".dat";
  std::ofstream file;
  file.open(filename);
  file << 0*h << "\t" << phi0 << std::endl; 
  if(file.is_open()){
    for(int k=0; k<N; ++k){
      phi = advanceRungeKutta(phi,method,h);
      file << (k+1)*h << "\t" << phi << std::endl;
    }
    std::cout << "Integration completed." << std::endl;
    return 1;
  }
  else
    std::cerr << "ERROR when opening file" << filename << std::endl;
  
  return 0;
}
