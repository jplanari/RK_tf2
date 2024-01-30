void predictor(tf2::Field &up, tf2::Field &un, tf2::Field &conv, tf2::Field &diff, tf2::Field &ms, double coef, double dt, tf2::Simulation &sim)
{
    tf2::oper_copy(un, up);
    tf2::oper_axpy(diff, up, dt*coef, 1.0);
    tf2::oper_axpy(conv, up, -dt*coef, 1.0);
    tf2::oper_axpy(ms, up, dt*coef, 1.0);
}

void predictorVector(tf2::Field &up, tf2::Field &un, std::vector<tf2::Field*> conv, std::vector<tf2::Field*> diff, tf2::Field &ms, butcherTableau coefs, double dt, int i, tf2::Simulation &sim)
{
  tf2::oper_copy(un,up);

  if(i!=0)
    for (int j=1; j<i; ++j)
      predictor(up,up,*(conv.at(j-1)),*(diff.at(j-1)),ms,coefs.A.at(id(i,j)),dt,sim);
  else
    for (int j=0; j<coefs.b.size(); ++j)
      predictor(up,up,*(conv.at(j)),*(diff.at(j)),ms,coefs.b.at(j),dt,sim);
} 

void projection(tf2::Field &uf, tf2::Field &p, tf2::Matrix &G)
{
    tf2::oper_prod(G,p,uf,-1.0,1.0);
}

void projection(tf2::Field &up, tf2::Field &p, tf2::Matrix &G, tf2::Matrix &IFC, tf2::Simulation &sim)
{
    auto &Gp = TF_getTmpField(sim,up.dim,"Faces");
    tf2::oper_prod(G,p,Gp);
    tf2::oper_prod(IFC,Gp,up,-1.0,1.0);
}

void pRHS(tf2::Field &pS, tf2::Field &ux, tf2::Field &uy, tf2::Field &uz, tf2::Simulation &sim)
{
    auto &Omega = tf2::getField(sim, "Omega_C");
    auto &DX    = tf2::getMatrix(sim, "DivX_FC");
    auto &DY    = tf2::getMatrix(sim, "DivY_FC");
    auto &DZ    = tf2::getMatrix(sim, "DivZ_FC");
    
    tf2::oper_prod(DX, ux, pS);
    tf2::oper_prod(DY, uy, pS, 1.0, 1.0);
    tf2::oper_prod(DZ, uz, pS, 1.0, 1.0);
    tf2::oper_prod(Omega, pS, -1.0);
}

void ID_CN_bocos(tf2::Field &phiC, tf2::Field &phiN, std::string fieldName, tf2::Simulation &sim)
{
  std::string bndsrcName = fieldName + "_BndSrc";
  std::string bndapplyName = fieldName + "_BndApply_NN";

  auto &aux = TF_getTmpField(sim, phiN);
  auto &phi_bnd = tf2::getField(sim, bndsrcName);
  auto &MB = tf2::getMatrix(sim, bndapplyName);
  auto &cm = tf2::meshNodeCellMask(sim);
  auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
  
  tf2::oper_prod(ID_CN, phiC, aux);
  tf2::oper_copy(phi_bnd, phiN);
  tf2::oper_fmadd(cm, aux, phiN, 1.0, 1.0);
  tf2::oper_prod(MB, aux, phiN, 1.0, 1.0);
}
