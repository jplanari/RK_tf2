#include <iostream>
#include <cmath>
#include <vector>
#include <string>
#include <fstream>
#include <complex>

#include "../../mods/fsm_rk/butcherTableaux.h"

std::complex<double> f(std::complex<double> phi, double angle)
{
  std::complex<double> lambda(-1.0*cos(angle*M_PI/180),sin(angle*M_PI/180));
  return lambda*phi;
}

std::complex<double> advanceRungeKutta(std::complex<double> phi, butcherTableau method, double h, double angle)
{
    std::complex<double> phin;
    int s = method.b.size();
    std::vector<std::complex<double>> phis(s);
    
    for(int i=0; i<s; ++i){
      phis.at(i) = phi;
      for(int j=0;j<i;++j)
        phis.at(i) += h*method.A.at(id(i+1,j+1))*f(phis.at(j),angle);
    }
    
    phin = phi;
    for(int i=0; i<s; ++i)
      phin += h*method.b.at(i)*f(phis.at(i),angle);

  return phin;
}

std::complex<double> advanceAdamsBashforth(std::complex<double> phi, std::complex<double> &phi0, double h, double angle)
{
  std::complex<double> phin;

  phin = phi+1.5*h*f(phi,angle)-0.5*h*f(phi0,angle);
  phi0=phi;

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
  double angle;
  
  if (caseName == "DIFFUSIVE") angle = 0.0;
  else if (caseName == "IMAGINARY") angle = 90.0;
  else if (caseName == "ANGLE") angle = atof(argv[5]);
  else {std::cerr << "Please introduce a caseName" << std::endl; return 1;}

  std::string filename = caseName+"/results/"+schemeName+"/testPositivity_"+schemeName+"_h_"+argv[2]+"_N_"+argv[3]+".dat";
  std::ofstream file;
  file.open(filename);
  file << 0*h << "\t" << phi0.real() << std::endl; 
  if(file.is_open()){
    for(int k=0; k<N; ++k){
      phi = advanceRungeKutta(phi,method,h,angle);
      file << (k+1)*h << "\t" << phi.real() << std::endl;
    }
    return 1;
  }
  else
    std::cerr << "ERROR when opening file" << filename << std::endl;
  
  return 0;
}
