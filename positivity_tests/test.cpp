#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <complex>

#include "../mods/fsm_rk/butcherTableaux.h"

std::complex<double> f(std::complex<double> phi, std::string casename)
{
  std::complex<double> lambda;
  if(casename == "IMAGINARY"){
    std::complex<double> lambdaI(0.0,1.0);
    lambda=lambdaI;
  }
  else if (casename == "DIFFUSIVE"){
    std::complex<double>lambdaD(-1.0,0.0);
    lambda=lambdaD;
  }
  
  return lambda*phi;
}

std::complex<double> advanceRungeKutta(std::complex<double> phi, butcherTableau method, double h, std::string casename)
{
    std::complex<double> phin;
    int s = method.b.size();
    std::vector<std::complex<double>> phis(s);
    
    for(int i=0; i<s; ++i){
      phis.at(i) = phi;
      for(int j=0;j<i;++j)
        phis.at(i) += h*method.A.at(id(i+1,j+1))*f(phis.at(j),casename);
    }
    
    phin = phi;
    for(int i=0; i<s; ++i)
      phin += h*method.b.at(i)*f(phis.at(i),casename);

  return phin;
}


int main(int argc, char** argv){
  std::string schemeName = argv[1];
  butcherTableau method = intScheme(schemeName);
//  displayMethod(method);
  double h = atof(argv[2]);
  std::complex<double> phi0(1.0,0.0);
  std::complex<double> phi=phi0;
  int N = atoi(argv[3]);
  std::string caseName = argv[4];
  std::string filename = caseName+"/results/"+schemeName+"/testPositivity_"+schemeName+"_h_"+argv[2]+"_N_"+argv[3]+".dat";
  std::ofstream file;
  file.open(filename);
  file << 0*h << "\t" << phi0 << std::endl; 
  if(file.is_open()){
    for(int k=0; k<N; ++k){
      phi = advanceRungeKutta(phi,method,h,caseName);
      file << (k+1)*h << "\t" << phi.real() << std::endl;
    }
    return 1;
  }
  else
    std::cerr << "ERROR when opening file" << filename << std::endl;
  
  return 0;
}
