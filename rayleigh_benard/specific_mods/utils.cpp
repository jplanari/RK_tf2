#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include <vector>

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
    double Ra = sim.IOParamD["Ra"]; //Rayleigh
    double Pr = sim.IOParamD["Pr"]; //Prandtl

    double &kinVisc = sim.IOParamD["kinVisc"];
    double &lambda = sim.IOParamD["lambda"];

    kinVisc = sqrt(Pr/Ra);
    lambda = 1/sqrt(Pr*Ra);
    
    sim.IOParamD["_MaxTime"] = 100;

    //setup_cfl(sim); 

    tf2::info("init_props completed.\n");
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

TF_Func bool monitor(tf2::Simulation &sim)
{
    static bool first = true;
    static bool steady_state = false;
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

