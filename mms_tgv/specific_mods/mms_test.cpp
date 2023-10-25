#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include <vector>
#include <cmath>
#include "../../mods/fsm_rk/butcherTableaux.h"

TF_Func void init_profile(tf2::Simulation &sim)
{
    auto dim = tf2::getField(sim, "ux_N").dim;
    auto profile = [=](double x, double y, double) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, sin(x)*cos(y));
        return result;
    };

    tf2::initField(sim, "ux_N", profile);

    auto profiley = [=](double x, double y, double) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, -1.0*cos(x)*sin(y));
        return result;
    };
    
    tf2::initField(sim, "uy_N", profiley);

    auto profilez = [=](double x, double y, double) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, -0.0*cos(x)*sin(y));
        return result;
    };
    
    tf2::initField(sim, "uz_N", profilez);
}

TF_Func void checkError(tf2::Simulation &sim)
{
  double t = sim.IOParamD["_ElapsedTime"];
  double Re = sim.IOParamD["Re"];

  auto &ux_anal = tf2::getOrCreateField(sim,1,"ux_N_anal","Nodes");
  auto &uy_anal = tf2::getOrCreateField(sim,1,"uy_N_anal","Nodes");
  
  auto profile = [=](double x, double y, double) -> tf2::Vn
  {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(ux_anal.dim, sin(x)*cos(y));
        return result;
  };

  tf2::initField(sim, "ux_N_anal", profile);
  
  tf2::oper_axpy(ux_anal,ux_anal,exp(-2*t/Re),0.0);

  auto profiley = [=](double x, double y, double) -> tf2::Vn
  {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(ux_anal.dim, -1.0*cos(x)*sin(y));
        return result;
  };
    
  tf2::initField(sim, "uy_N_anal", profiley);
  tf2::oper_axpy(uy_anal,uy_anal,exp(-2*t/Re),0.0);
  
  auto &ux = tf2::getField(sim,"ux_N");
  auto &uy = tf2::getField(sim,"uy_N");

  tf2::oper_axpy(ux,ux_anal,-1.0,1.0);
  tf2::oper_apply(ux_anal, [](double x) {return fabs(x);});

  tf2::oper_axpy(uy,uy_anal,-1.0,1.0);
  tf2::oper_apply(uy_anal, [](double x) {return fabs(x);});
  
  tf2::info("maxU=%e, maxV=%e\n",tf2::oper_max(ux_anal),tf2::oper_max(uy_anal));
  tf2::info("t=%e\tERROR Ux = %e\tERROR Uy = %e\n",t,tf2::oper_max(ux_anal)[0],tf2::oper_max(uy_anal)[0]);
}

TF_Func bool mmsRK(tf2::Simulation &sim)
{
  if(sim.IOParamD["_ElapsedTime"] > sim.IOParamD["_MaxTime"]) return tf2::Iter_Stop; 
  
  auto &ux_N = tf2::getField(sim,"ux_N");
  auto &uy_N = tf2::getField(sim,"uy_N");

  auto &uxp = TF_getTmpField(sim,ux_N);
  auto &uyp = TF_getTmpField(sim,uy_N);

  double Re = sim.IOParamD["Re"];
  
  butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

  std::vector<double> b = coefs.b;
  std::vector<double> A = coefs.A;
  int s = b.size();

  std::vector<tf2::Field*> ux(s);
  std::vector<tf2::Field*> uy(s);
   
  ux.at(0) = &uxp;
  uy.at(0) = &uyp;

  tf2::oper_axpy(ux_N,uxp,-2/Re,0.0);
  tf2::oper_axpy(uy_N,uyp,-2/Re,0.0);
 
  for(int i=1; i < s; ++i)
  {
    auto &uxi = tf2::getOrCreateField(sim,ux_N.dim,"ux_"+std::to_string(i),"Nodes");
    auto &uyi = tf2::getOrCreateField(sim,ux_N.dim,"uy_"+std::to_string(i),"Nodes");
    
    tf2::oper_copy(ux_N, uxi);
    tf2::oper_copy(uy_N, uyi);

    for(int j=1; j<i; ++j)
    {
      tf2::oper_axpy(*(ux.at(j-1)),uxi,A.at(id(i,j)),1.0);
      tf2::oper_axpy(*(uy.at(j-1)),uyi,A.at(id(i,j)),1.0); 
    }
    
    tf2::oper_axpy(uxi,uxi,-2/Re,0.0);
    tf2::oper_axpy(uyi,uyi,-2/Re,0.0);

    ux.at(i) = &uxi;
    uy.at(i) = &uyi; 
  } 

  
  for(int i=0; i < s; ++i)
  {
      tf2::oper_axpy(*(ux.at(i)),ux_N,b.at(i),1.0);
      tf2::oper_axpy(*(uy.at(i)),uy_N,b.at(i),1.0); 
  }

  checkError(sim);

  return tf2::Iter_Continue;

}


