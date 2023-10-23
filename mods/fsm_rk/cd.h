
void diffusive(tf2::Field &u, tf2::Field &diff, tf2::Simulation &sim)
{
    auto &L = tf2::getMatrix(sim, "Lap_NC");

    TF_uAssert(sim.IOParamD.count("kinVisc") > 0, "Parameter 'kinVisc' not defined");
    double visc = sim.IOParamD["kinVisc"];

    // Diffusive_i = visc*lap(u_i)
    tf2::oper_prod(L, u, diff, visc);
}

void convective(tf2::Field &u, tf2::Field &conv, tf2::Simulation &sim)
{
    auto &DX  = tf2::getMatrix(sim, "DivX_FC");
    auto &DY  = tf2::getMatrix(sim, "DivY_FC");
    auto &DZ  = tf2::getMatrix(sim, "DivZ_FC");

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
}


