
tf2::Field &diffusive(tf2::Field &u, uint32_t i, tf2::Simulation &sim)
{
    auto &L = tf2::getMatrix(sim, "Lap_NC");

    std::string name = "diff_" + std::to_string(i);

    auto &diff = tf2::getOrCreateField(sim,u.dim,name,"Cells");

    TF_uAssert(sim.IOParamD.count("kinVisc") > 0, "Parameter 'kinVisc' not defined");
    double visc = sim.IOParamD["kinVisc"];

    // Diffusive_i = visc*lap(u_i)
    tf2::oper_prod(L, u, diff, visc);

    return diff;
}

tf2::Field &convective(tf2::Field &u, uint32_t i, tf2::Simulation &sim)
{
    auto &DX  = tf2::getMatrix(sim, "DivX_FC");
    auto &DY  = tf2::getMatrix(sim, "DivY_FC");
    auto &DZ  = tf2::getMatrix(sim, "DivZ_FC");

    std::string name = "conv_" + std::to_string(i);

    auto &conv = tf2::getOrCreateField(sim,u.dim,name,"Cells");

    auto &ufx   = tf2::getField(sim, "ux_F");
    auto &ufy   = tf2::getField(sim, "uy_F");
    auto &ufz   = tf2::getField(sim, "uz_F");

    // Temporary variables.
    auto &aux = TF_getTmpField(sim, ufx);
    // Convective_i = div(u*u_i)

    
    tf2::oper_axpy(ufx,aux,1.0,0.0);
    tf2::oper_prod(u,aux,1.0);
    tf2::oper_prod(DX, aux, conv);
    tf2::oper_axpy(ufy,aux,1.0,0.0);
    tf2::oper_prod(u,aux,1.0);
    tf2::oper_prod(DY, aux, conv, 1.0, 1.0);
    tf2::oper_axpy(ufz,aux,1.0,0.0);
    tf2::oper_prod(u,aux,1.0);
    tf2::oper_prod(DZ, aux, conv, 1.0, 1.0);
    
    return conv;
}


