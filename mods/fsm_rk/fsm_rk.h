tf2::Field &predictor(tf2::Field &un, tf2::Field &conv, tf2::Field &diff, tf2::Field &ms, double coef, double dt, tf2::Simulation &sim)
{
    auto &up = TF_getTmpField(sim,un);

    tf2::oper_copy(un, up);
    tf2::oper_axpy(diff, up, dt*coef, 1.0);
    tf2::oper_axpy(conv, up, -dt*coef, 1.0);
    tf2::oper_axpy(ms, up, dt*coef, 1.0);

    return up;
}

tf2::Field &predictorVector(tf2::Field &un, std::vector<tf2::Field*> conv, std::vector<tf2::Field*> diff, tf2::Field &ms, butcherTableau coefs, double dt, int i, tf2::Simulation &sim)
{
  auto &up = TF_getTmpField(sim,un);
  tf2::oper_copy(un,up);

  if(i!=0)
    for (int j=1; j<i; ++j)
      up = predictor(up,*(conv.at(j-1)),*(diff.at(j-1)),ms,coefs.A.at(id(i,j)),dt,sim);
  else
    for (int j=0; j<coefs.b.size(); ++j)
      up = predictor(up,*(conv.at(j)),*(diff.at(j)),ms,coefs.b.at(j),dt,sim);
  
  return up;
} 

void projection(tf2::Field &uf, tf2::Field &p, tf2::Matrix &G)
{
    tf2::oper_prod(G,p,uf,-1.0,1.0);
}

void projection(tf2::Field &up, tf2::Field &p, tf2::Matrix &G, tf2::Matrix &IFC, tf2::Simulation &sim)
{
    auto &Gp = TF_getTmpField(sim,1,"Faces");
    tf2::oper_prod(G,p,Gp);
    tf2::oper_prod(IFC,Gp,up,-1.0,1.0);
}

tf2::Field &pRHS(tf2::Field &ux, tf2::Field &uy, tf2::Field &uz, tf2::Simulation &sim)
{
    auto &pS    = tf2::getOrCreateField(sim,ux.dim,"pSource","Cells");
    auto &Omega = tf2::getField(sim, "Omega_C");
    auto &DX    = tf2::getMatrix(sim, "DivX_FC");
    auto &DY    = tf2::getMatrix(sim, "DivY_FC");
    auto &DZ    = tf2::getMatrix(sim, "DivZ_FC");
    
    tf2::oper_prod(DX, ux, pS);
    tf2::oper_prod(DY, uy, pS, 1.0, 1.0);
    tf2::oper_prod(DZ, uz, pS, 1.0, 1.0);
    tf2::oper_prod(Omega, pS, -1.0);

    return pS;
}


