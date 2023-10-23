#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include <vector>

TF_Func void init_props(tf2::Simulation &sim)
{
    double Re_tau = sim.IOParamD["Re"];
    double &kinVisc = sim.IOParamD["kinVisc"];
    double Lx = tf2::getLength(tf2::getSMesh(sim))[0];
    double U0 = sim.IOParamD["U0"];

    sim.IOParamD["_MaxTime"] = (Lx/U0)*20.0;

    kinVisc = 1.0/Re_tau;
}

double my_rand0(double)
{
    return (1.0-2.0*drand48());
}

TF_Func void init_profile(tf2::Simulation &sim)
{
    auto U0 = sim.IOParamD["U0"];
    auto dim = tf2::getField(sim, "ux_N").dim;
    
    auto profile_x = [=](double x, double y, double z) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, U0*sin(x)*cos(y)*cos(z));
        return result;
    };
 
    auto profile_y = [=](double x, double y, double z) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, U0*cos(x)*sin(y)*cos(z));
        return result;
    };   

    auto profile_z = [=](double x, double y, double z) -> tf2::Vn
    {
        // We want to apply the same profile for all components / simulations,
        // so we just replicate it 'dim' times.
        tf2::Vn result(dim, U0*cos(x)*cos(y)*sin(z));
        return result;
    };   

    tf2::initField(sim, "ux_N", profile_x);
    tf2::initField(sim, "uy_N", profile_y);
    tf2::initField(sim, "uz_N", profile_z);

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

