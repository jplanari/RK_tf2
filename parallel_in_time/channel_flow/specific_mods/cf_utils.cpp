#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include <vector>
#include <fstream>

void setup_cfl(tf2::Simulation &sim)
{
    double &DX = sim.IOParamD["dx"];
    auto &omega = tf2::getField(sim, "Omega_C");
    double v = tf2::oper_min(omega)[0];
    DX = tf2::allReduce(pow(v, 1.0/3.0), MPI_MIN);
    const double kinVisc = sim.IOParamD["kinVisc"];
    double &tVisc = sim.IOParamD["tVisc"];
    tVisc = 0.2*DX*DX/kinVisc;
    sim.IOParamD["_TimeStep"] = tVisc;
}

TF_Func void init_props(tf2::Simulation &sim)
{
    double Re_tau = sim.IOParamD["Re_tau"];
    double &kinVisc = sim.IOParamD["kinVisc"];
    double &Ub = sim.IOParamD["Ub"];
    double &d = sim.IOParamD["delta"];
    double Ly = tf2::getLength(tf2::getSMesh(sim))[1];
    sim.IOParamD["_MaxTime"] = 1e5;
    sim.IOParamD["_RAn"] = Ub;
    sim.IOParamD["_RAn1"] = Ub;

    d = 0.5*Ly;
    Ub = d*0.5*Re_tau;
    kinVisc = 1.0/Re_tau;
    auto &msx = tf2::getField(sim, "momSrcx_C");
    double gp = sim.IOParamD["gradp"];
    tf2::oper_setConst(msx, gp);

    setup_cfl(sim);
    tf2::info("init_props completed.\n");
}

double my_rand0(double)
{
    return (1.0-2.0*drand48());
}

TF_Func void randomize_u(tf2::Simulation &sim)
{
    srand48(tf2::mpiRank());
    auto &uy = tf2::getField(sim, "uy_N");
    auto &uz = tf2::getField(sim, "uz_N");

    // Since we are passing 'my_rand0' as the argument to 'oper_apply', and
    // 'my_rand0' is a function that maps a double into a double (as opposed to
    // mapping a tf2::Vn into a tf2::Vn), it will work irrespective of how many
    // components the field has: it will be applied for every component of
    // every element.
    tf2::oper_apply(uy, my_rand0);
    tf2::oper_apply(uz, my_rand0);
    tf2::info("randomize_u completed\n");
}

TF_Func void init_profile_ux(tf2::Simulation &sim)
{
    auto Ub = sim.IOParamD["Ub"];
    auto d = sim.IOParamD["delta"];
    auto dim = tf2::getField(sim, "ux_N").dim;
    auto profile = [=](double, double y, double) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, Ub*y/d*(2.0-y/d)+0.01*Ub*my_rand0(.1));
        return result;
    };
    tf2::initField(sim, "ux_N", profile);
    tf2::info("init_profile_ux completed.\n");

    auto &INF = tf2::getMatrix(sim,"Interp_NF");

    auto &ux = tf2::getField(sim,"ux_N");
    auto &uy = tf2::getField(sim,"uy_N");
    auto &uz = tf2::getField(sim,"uz_N");

    auto &uxf = tf2::getField(sim,"ux_F");
    auto &uyf = tf2::getField(sim,"uy_F");
    auto &uzf = tf2::getField(sim,"uz_F");

    tf2::oper_prod(INF,ux,uxf);
    tf2::oper_prod(INF,uy,uyf);
    tf2::oper_prod(INF,uz,uzf);
}

void printPhi(tf2::Simulation &sim)
{
  std::string Re = std::to_string(sim.IOParamD["Re_tau"]);
  std::string Rkfct = std::to_string(sim.IOParamD["RKfct"]);
  std::string RKmethod = sim.IOParamS["RKmethod"];

  std::string filename = "results/phi_Re_tau-"+Re+"_RKfct-"+Rkfct+"_"+RKmethod+".dat";
  std::ofstream file;

  static bool first = true;
  if (first) {
    file.open(filename);
    first = false;
  }
  else
    file.open(filename,std::ios::app);
  
  if (file.is_open()){
    file << sim.IOParamI["_Iter"] << "\t" << atan(sim.IOParamD["_EVimag"]/sim.IOParamD["_EVreal"]) << std::endl; 
    file.close();;
  }
}

TF_Func void init_omega(tf2::Simulation &sim)
{
    // For this particular application we need a field with the volume of each
    // cell, replicated in all the components. This volume field will be used
    // in the FSM, to scale the divergence of the predictor velocity. This is
    // needed because we use the laplacian multiplied by the volume as the
    // pressure equation, in order to make the matrix SPD.
    auto dim    = tf2::getField(sim, "ux_N").dim;
    auto &omega = tf2::getOrCreateField(sim, dim, "Omega_C", "Cells");
    const auto &m = tf2::getSMesh(sim);
    auto cells = tf2::getNumCells(m);
    std::vector<double> buffer(tf2::getNumEntries(omega));

    for (uint32_t c = 0; c < cells; c++)
    {
        const auto &C = tf2::getCellIJKFromId(m, c);
        double v = tf2::calcCellVolume(m, C);
        for (uint32_t d = 0; d < dim; d++)
        {
            buffer[c*dim+d] = v;
        }
    }
    tf2::oper_setData(omega, buffer.data());
    tf2::info("init_omega completed.\n");
}

TF_Func bool update_DT_cfl(tf2::Simulation &sim) 
{
	auto &ux  = tf2::getField(sim, "ux_N");
	auto &uy  = tf2::getField(sim, "uy_N");
	auto &uz  = tf2::getField(sim, "uz_N");
	auto &mod = TF_getTmpField(sim, ux);
	auto &tmp = TF_getTmpField(sim, ux);

	tf2::oper_axpy(ux,mod,1.0,0.0);
	tf2::oper_prod(ux,mod,1.0);
	
	tf2::oper_axpy(uy,tmp,1.0,0.0);
	tf2::oper_prod(uy,tmp,1.0);
	tf2::oper_axpy(tmp,mod,1.0,1.0);

	tf2::oper_axpy(uz,tmp,1.0,0.0);
	tf2::oper_prod(uz,tmp,1.0);
	tf2::oper_axpy(tmp,mod,1.0,1.0);
	
	const double uref = sqrt(tf2::max(tf2::oper_max(mod)));
	const double cfl = sim.IOParamD["cfl"];
	const double tVisc = sim.IOParamD["tVisc"];
	const double DX = sim.IOParamD["dx"]; 
	sim.IOParamD["_TimeStep"] = cfl*std::min(0.35*DX/uref, tVisc);
	return tf2::Iter_Continue;
} 

bool determineSteadyState(double max_x, double max_u, tf2::Simulation &sim)
{
  if(sim.IOParamI["_Iter"]%1000 == 0)
  {
      if(fabs(max_x-max_u)/max_u < 1e-1){
        sim.IOParamD["_MaxTime"] = sim.IOParamD["_ElapsedTime"] + sim.IOParamD["nFT"]*4*M_PI/(0.65*sim.IOParamD["maxU"]);
        sim.IOParamI["_TAVG_Start"] = sim.IOParamI["_Iter"];
        return true;
      }
      sim.IOParamD["maxU"] = max_x;
      return false;
  } 
}

bool rollingSteadyState(double max_x, int window_size, tf2::Simulation &sim)
{
  double curr_ra = sim.IOParamD["_RAn"];
  curr_ra += max_x;
  if(sim.IOParamI["_Iter"]%window_size==0)
  {
    curr_ra /= (double) window_size;
    if (fabs(curr_ra-sim.IOParamD["_RAn1"])/sim.IOParamD["_RAn1"] < 1e-1){
        sim.IOParamD["_MaxTime"] = sim.IOParamD["_ElapsedTime"] + sim.IOParamD["nFT"]*4*M_PI/(0.65*sim.IOParamD["maxU"]);
        sim.IOParamI["_TAVG_Start"] = sim.IOParamI["_Iter"];
        return true;
    }
    else
      sim.IOParamD["_RAn1"] = sim.IOParamD["_RAn"];
  }
  return false;
}

TF_Func bool monitor(tf2::Simulation &sim)
{
    static bool first = true;
    const uint32_t ndim = tf2::getField(sim,"ux_N").dim;
    if (ndim>1){
    if (first)
    {
        tf2::info("# "
                  " ite wclock "
                  " ts time "
                  " resSolver "
                  "\n");
        first = false;
    }

    if (sim.IOParamI["_Iter"]%10 == 0)
    {
    static auto s = tf2::getSolver(sim, "Pressure_Solver");
    tf2::info("real=%e, imag=%e\n",sim.IOParamD["_EVreal"],sim.IOParamD["_EVimag"]);
    tf2::info("%d %.8e %.5e %.5e %.5e\n",
              sim.IOParamI["_Iter"], runTime(sim),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tf2::oper_solve_residual(s));
    }
    }
    else
    {
    if (first)
    {
        tf2::info("# "
                  " ite wclock "
                  " ts time "
                  " resSolver "
		  " max(ux) "
                  "\n");
        first = false;
    }
    double max_x;
    if (sim.IOParamI["_Iter"]%10 == 0)
    {
    printPhi(sim);
    auto &ux = tf2::getField(sim, "ux_N");
    max_x = std::max(fabs(tf2::oper_max(ux)[0]), fabs(tf2::oper_min(ux)[0]));
    static auto s = tf2::getSolver(sim, "Pressure_Solver");
    tf2::info("%d %.8e %.5e %.5e %.5e %.5e\n",
              sim.IOParamI["_Iter"], runTime(sim),
              sim.IOParamD["_TimeStep"], sim.IOParamD["_ElapsedTime"],
              tf2::oper_solve_residual(s),max_x);
    }
    }
    return tf2::Iter_Continue;
}

TF_Func void setupManualBocos(tf2::Simulation &sim)
{
auto &u = tf2::getField(sim, "ux_N");
auto cells = tf2::getDomainSize(sim, "Cells");
std::vector<double> buffer(tf2::getNumEntries(u), 1.0);
for (uint32_t it = cells*u.dim; it < buffer.size(); it++)
{
    buffer[it] = 0.0;
}
auto &cellMask = tf2::getOrCreateField(sim, "bnd", u);
tf2::oper_setData(cellMask, buffer.data());
}

TF_Func bool applyManualBocos(tf2::Simulation &sim)
{
    auto &bnd = tf2::getField(sim,"bnd");
	
    auto &ux = tf2::getField(sim,"ux_N");
    auto &uy = tf2::getField(sim,"uy_N");
    auto &uz = tf2::getField(sim,"uz_N");

    tf2::oper_prod(bnd,ux);
    tf2::oper_prod(bnd,uy);
    tf2::oper_prod(bnd,uz);

    return tf2::Iter_Continue;
}

TF_Func void calc_div_u(tf2::Simulation &sim)
{
    auto &DX = tf2::getMatrix(sim, "DivX_FC");
    auto &DY = tf2::getMatrix(sim, "DivY_FC");
    auto &DZ = tf2::getMatrix(sim, "DivZ_FC");

    auto &ux = tf2::getField(sim, "ux_F");
    auto &uy = tf2::getField(sim, "uy_F");
    auto &uz = tf2::getField(sim, "uz_F");

    auto &div = tf2::getOrCreateField(sim, "div(u)", DX, ux);
    tf2::oper_prod(DX, ux, div, 1.0, 0.0);
    tf2::oper_prod(DY, uy, div, 1.0, 1.0);
    tf2::oper_prod(DZ, uz, div, 1.0, 1.0);
}

