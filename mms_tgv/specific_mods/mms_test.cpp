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
        tf2::Vn result(dim, 0.5);
        return result;
    };

    tf2::initField(sim, "ux_N", profile);

    auto profiley = [=](double x, double y, double) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, -0.5);
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
        tf2::Vn result(ux_anal.dim, 0.5);
        return result;
  };

  tf2::initField(sim, "ux_N_anal", profile);
  
  tf2::oper_axpy(ux_anal,ux_anal,exp(-2*t/Re),0.0);

  auto profiley = [=](double x, double y, double) -> tf2::Vn
  {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(ux_anal.dim, -0.5);
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

  tf2::info("ux=%.10f\n",tf2::oper_max(ux)[0]);
  tf2::info("t=%e\tERROR Ux = %e\nERROR Uy = %e\n",t,tf2::oper_max(ux_anal)[0],tf2::oper_max(uy_anal)[0]);
}

TF_Func bool mmsRK(tf2::Simulation &sim)
{
  if(sim.IOParamD["_ElapsedTime"] > sim.IOParamD["_MaxTime"]) return tf2::Iter_Stop; 
  
  auto &ux_N = tf2::getField(sim,"ux_N");
  auto &uy_N = tf2::getField(sim,"uy_N");

  double Re = sim.IOParamD["Re"];

  double dt = sim.IOParamD["_TimeStep"];

  double f = -2.0/Re;
  
  butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

  std::vector<double> b = coefs.b;
  std::vector<double> A = coefs.A;
  int s = b.size();

  std::vector<tf2::Field*> ux(s);
  std::vector<tf2::Field*> uy(s);
   
  ux.at(0) = &ux_N;
  uy.at(0) = &uy_N;

  for(int i=2; i <= s; ++i)
  {
    auto &uxi = tf2::getOrCreateField(sim,ux_N.dim,"ux_"+std::to_string(i),"Nodes");
    auto &uyi = tf2::getOrCreateField(sim,ux_N.dim,"uy_"+std::to_string(i),"Nodes");
    
    tf2::oper_copy(ux_N, uxi);
    tf2::oper_copy(uy_N, uyi);

    for(int j=1; j<i; ++j) //sum for j=1 to i-1
    {
      tf2::oper_axpy(*(ux.at(j-1)),uxi,dt*f*A.at(id(i,j)),1.0);
      tf2::oper_axpy(*(uy.at(j-1)),uyi,dt*f*A.at(id(i,j)),1.0); 
    }
  
    ux.at(i-1) = &uxi;
    uy.at(i-1) = &uyi; 
  
  } 
  
  for(int i=0; i < s; ++i)
  {
      tf2::oper_axpy(*(ux.at(i)),ux_N,dt*f*b.at(i),1.0);
      tf2::oper_axpy(*(uy.at(i)),uy_N,dt*f*b.at(i),1.0); 
  }

  return tf2::Iter_Continue;

}


