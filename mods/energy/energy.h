void generateSourceTerm(double alpha, tf2::Field &T, tf2::Field &source, tf2::Simulation &sim)
{
  auto &Tc = TF_getTmpField(sim,source);
  auto &M = tf2::getMatrix(sim, "Id_NC");
  tf2::oper_prod(M,T,Tc);
  tf2::oper_axpy(Tc,source,alpha,0.0);
}


