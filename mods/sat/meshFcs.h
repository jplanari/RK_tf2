tf2::Field &calcFaceStaggeredVolumes(tf2::Simulation &sim)
{
    if(tf2::hasField(sim, "_FaceVolume")) return tf2::getField(sim, "_FaceVolume");
    
    auto dim = tf2::getField(sim,"ux_N").dim;
    auto &result = tf2::newField(sim, dim, "_FaceVolume", "Faces");
    
    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);

    if (tf2::hasSMesh(sim))
    {
      const auto &m = tf2::getSMesh(sim);
      for (uint32_t it = 0; it < resultSize/dim; it++)
      {
        const auto &F = tf2::getFaceIJKFromId(m,it);
        for(int d=0; d<dim; d++)
        buffer[it*dim+d] = tf2::calcFaceStaggeredVol(m,F);
      }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

tf2::Field &faceCellDist(tf2::Simulation &sim)
{
  if(tf2::hasField(sim,"_FaceCellDist")) return tf2::getField(sim,"_FaceCellDist");

  auto dim = tf2::getField(sim,"ux_N").dim;
  auto &result = tf2::newField(sim, dim, "_FaceCellDist", "Faces");
  uint32_t resultSize = tf2::getNumEntries(result);
  std::vector<double> buffer(resultSize);

  if(tf2::hasSMesh(sim))
  {
    const auto &m = tf2::getSMesh(sim);
    for (uint32_t it = 0; it < resultSize/dim; it++)
    {
      const auto F = tf2::getFaceIJKFromId(m, it);
      const auto n = tf2::calcFaceNormal(m, F);
      if (tf2::isBoundaryFace(m, F))
      {
        const auto A = tf2::calcCellDim(m, F);
        for(int d=0; d<dim; d++)
            buffer[it*dim+d] = 0.5*fabs(A*n);
      }
      else
      {
        const auto O = tf2::getOtherCellIJKFromFace(m, F);
        for(int d=0; d<dim; d++)
            buffer[it*dim+d] = tf2::calcCellDist(m, F.i, F.j, F.k, O.i, O.j, O.k);
      }
    } 
  }

  tf2::oper_setData(result, buffer.data());
  return result;
}



tf2::IJK returnNbIJK(tf2::IJK CNB, uint32_t f)
{
  tf2::IJK C = CNB;
  switch(f){
    case 0:
      C.i--;
      break;
    case 1:
      C.i++;
      break;
    case 2:
      C.j--;
      break;
    case 3:
      C.j++;
      break;
    case 4:
      C.k--;
      break;
    case 5:
      C.k++;
      break;
  }
  return C;
}

std::vector<tf2::IJK> getFaceNbCells(tf2::IJK F, const tf2::SMesh &m, tf2::Simulation &sim)
{
  std::vector<tf2::IJK> result;
  result.push_back(F);
  uint32_t f = F.f;
  //If face F is not a boundary face:
  if(!tf2::isBoundaryFace(m,F)){
     result.push_back(returnNbIJK(F,f));}
  else
  {
    if(f%2!=0){ f--; result.at(0) = returnNbIJK(F,f);}
    else result.at(0) = F;
  }

  return result; 
}
