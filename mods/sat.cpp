#include "tf2/Simulation.h"
#include "tf2/Opers.h"
#include <cmath>

#include "fsm_rk/butcherTableaux.h"
#include "sat/stability.h"
#include "sat/eigenbounds.h"

void computeDT(tf2::Simulation &sim)
{
    double ev_r = sim.IOParamD["_EVreal"];
    double ev_i = sim.IOParamD["_EVimag"];
    double ev_norm = sqrt(ev_r*ev_r+ev_i*ev_i);

    int ns_rk = sim.IOParamI["order_RK"];
    double phi = M_PI-atan(ev_i/ev_r);

    sim.IOParamD["_TimeStep"] = sim.IOParamD["RKfct"]*stabilityRegion(phi,intScheme(sim.IOParamS["RKmethod"]))/ev_norm;
}

void computeDT_efficiency(tf2::Simulation &sim)
{
    double ev_r = sim.IOParamD["_EVreal"];
    double ev_i = sim.IOParamD["_EVimag"];
    double ev_norm = sqrt(ev_r*ev_r+ev_i*ev_i);
    
    if(sim.IOParamS["RKmethod"] != "efficiency") computeDT(sim);
    else{
      int idd;
      double phi = M_PI-atan(ev_i/ev_r);
      double h;
      double hmax=0.0, htilde, hh;
      butcherTableau candidate;
      double f;
    
      if(sim.IOParamS["FPJ"] == "yes") f = 0.85;
      else f = 0.0;

      std::vector<std::string> pool = multiRK_method(sim);
	    
      for(int id=0; id<pool.size(); ++id){
        candidate = intScheme(pool.at(id));
	      h = stabilityRegion(phi,candidate);
	      htilde = (h/ev_norm)/tauMethod(candidate.b.size(),f);
	      if(htilde>hmax){
		      hmax = htilde;
          hh = h/ev_norm;
		      idd = id;
	      }
	    }
      
      tf2::info("Scheme selected = %d\n", idd);
      sim.IOParamI["_schEff"] = idd; //idea is to remove this option, kept now just to check if everything works fine
      sim.IOParamS["RKname"] = pool.at(idd);
      sim.IOParamD["_TimeStep"] = sim.IOParamD["RKfct"]*hh;

    }
}

TF_Func void SetUp_SAT_Gershgorin_efficiency(tf2::Simulation &sim)
{
  computeEV_Gershgorin(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_Gershgorin_efficiency(tf2::Simulation &sim)
{
	computeImagEV_Gershgorin(sim);
	computeDT_efficiency(sim);
	return tf2::Iter_Continue;
}

TF_Func void SetUp_SAT_GershgorinMat_efficiency(tf2::Simulation &sim)
{
  computeEV_GershgorinMat(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_GershgorinMat_efficiency(tf2::Simulation &sim)
{
	computeImagEV_GershgorinMat(sim);
	computeDT_efficiency(sim);
	return tf2::Iter_Continue;
}
TF_Func void SetUp_SAT_AlgEigCD_efficiency(tf2::Simulation &sim)
{
  computeEV_AlgEigCD(sim);
	computeDT_efficiency(sim);
}

TF_Func bool Iter_SAT_AlgEigCD_efficiency(tf2::Simulation &sim)
{
	computeImagEV_AlgEigCD(sim);
	computeDT_efficiency(sim);
	return tf2::Iter_Continue;
}
