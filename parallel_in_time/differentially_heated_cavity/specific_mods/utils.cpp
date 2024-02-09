#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include <vector>

TF_Func void init_props(tf2::Simulation &sim)
{
    double Ra = sim.IOParamD["Ra"]; //Rayleigh
    double Pr = sim.IOParamD["Pr"]; //Prandtl

    double &kinVisc = sim.IOParamD["kinVisc"];
    double &lambda = sim.IOParamD["lambda"];

    kinVisc = sqrt(Pr/Ra);
    lambda = 1/sqrt(Pr*Ra);
    

    tf2::info("TAVG_start=%d\n",sim.IOParamI["_TAVG_Start"]);
    
    tf2::info("pre: max iters=%d\n",sim.cfgRT.maxIters);    
    sim.cfgRT.maxIters = sim.IOParamI["_TAVG_Start"] + (sim.cfgRT.maxIters - sim.IOParamI["_TAVG_Start"])/sim.IOParamI["nSims"];
    tf2::info("post: max iters=%d\n",sim.cfgRT.maxIters);    

    sim.IOParamD["_MaxTime"] = 100;

    if(sim.IOParamI["nSims"] > 1)
      tf2::info("Parallel-in-time, with %d rhs. Maximum number of iterations: %d\n",sim.IOParamI["nSims"],sim.cfgRT.maxIters);

    tf2::info("init_props completed.\n");
}

double my_rand0(double)
{
    return (1e-1-2e-1*drand48());
}

TF_Func void init_fields(tf2::Simulation &sim)
{
    srand48(tf2::mpiRank());
    auto &ux = tf2::getField(sim, "ux_N");
    auto &uy = tf2::getField(sim, "uy_N");
    auto &uz = tf2::getField(sim, "uz_N");

    // Since we are passing 'my_rand0' as the argument to 'oper_apply', and
    // 'my_rand0' is a function that maps a double into a double (as opposed to
    // mapping a tf2::Vn into a tf2::Vn), it will work irrespective of how many
    // components the field has: it will be applied for every component of
    // every element.
    tf2::oper_apply(ux, my_rand0);
    tf2::oper_apply(uy, my_rand0);
    tf2::oper_apply(uz, my_rand0);

    auto &INF = tf2::getMatrix(sim, "Interp_NF");
    auto &ufx = tf2::getField(sim, "ux_F");
    auto &ufy = tf2::getField(sim, "uy_F");
    auto &ufz = tf2::getField(sim, "uz_F");

    tf2::oper_prod(INF,ux,ufx);
    tf2::oper_prod(INF,uy,ufy);
    tf2::oper_prod(INF,uz,ufz);

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

