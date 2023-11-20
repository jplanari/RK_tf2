double calcFaceArea(const tf2::SMesh &m, tf2::IJK F)
{
    double result;
    auto areas = tf2::calcCellFacesArea(m, F);
    auto n = tf2::calcFaceNormal(m, F);
    result = fabs(areas*n);
    return result;
}

tf2::Field &getFaceNc(tf2::Simulation &sim, uint32_t comp, const std::string &meshName = "")
{
    const std::string name = "_FaceNormal"+std::to_string(comp);
    if (tf2::hasField(sim, name)) return tf2::getField(sim, name);
    const auto dim = tf2::getField(sim,"ux_N").dim;
    const auto &msuf  = tf2::meshSuffix(meshName);
    const auto &domain = "Faces"+msuf;
    auto &result = tf2::newField(sim, dim, name, domain);

    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);
    if (tf2::hasSMesh(sim))
    {
        const auto &m = tf2::getSMesh(sim);
        for (uint32_t it = 0; it < resultSize; it++)
        {
            const auto &F = tf2::getFaceIJKFromId(m, it);
            buffer[it] = tf2::calcFaceNormal(m, F)[comp];
        }
    }
    else
    {
        const auto &m = tf2::getUMesh(sim, meshName);
        uint32_t it = 0;
        for (auto f : m.getFacesRange())
        {
            buffer[it++] = m.getFaceNormal(f)[comp];
        }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

tf2::Field &getFaceNx(tf2::Simulation &sim, const std::string &meshName = "")
{
    return getFaceNc(sim, 0, meshName);
}

tf2::Field &getFaceNy(tf2::Simulation &sim, const std::string &meshName = "")
{
    return getFaceNc(sim, 1, meshName);
}

tf2::Field &getFaceNz(tf2::Simulation &sim, const std::string &meshName = "")
{
    return getFaceNc(sim, 2, meshName);
}

tf2::Field &getFaceArea(tf2::Simulation &sim, const std::string &meshName = "")
{
    const std::string &name = "_FaceArea";
    if (tf2::hasField(sim, name)) return tf2::getField(sim, name);

    const auto &msuf  = tf2::meshSuffix(meshName);
    const auto &domain = "Faces"+msuf;
    const auto dim = tf2::getField(sim,"ux_N").dim;

    auto &result = tf2::newField(sim, dim, name, domain);

    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);
    if (tf2::hasSMesh(sim))
    {
        const auto &m = tf2::getSMesh(sim);
        for (uint32_t it = 0; it < resultSize; it++)
        {
            const auto &F = tf2::getFaceIJKFromId(m, it);
            buffer[it] = calcFaceArea(m, F);
        }
    }
    else
    {
        const auto &m = tf2::getUMesh(sim, meshName);
        uint32_t it = 0;
        for (auto f : m.getFacesRange())
        {
            buffer[it++] = m.getFaceArea(f);
        }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

tf2::Field &getFaceInnerMask(tf2::Simulation &sim, const std::string &meshName = "")
{
    if (tf2::hasField(sim, "_FaceInnerMask")) return tf2::getField(sim, "_FaceInnerMask");

    const auto &msuf  = tf2::meshSuffix(meshName);
    const auto &domain = "Faces"+msuf;
    const auto dim = tf2::getField(sim,"ux_N").dim;
    auto &result = tf2::newField(sim, dim, "_FaceInnerMask", domain);

    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);
    if (tf2::hasSMesh(sim))
    {
        const auto &m = tf2::getSMesh(sim);
        for (uint32_t it = 0; it < resultSize; it++)
        {
            const auto &F = tf2::getFaceIJKFromId(m, it);
            buffer[it] = TO_F64(not tf2::isBoundaryFace(m, F));
        }
    }
    else
    {
        const auto &m = tf2::getUMesh(sim, meshName);
        uint32_t it = 0;
        for (auto f : m.getFacesRange())
        {
            buffer[it++] = TO_F64(not m.isBoundaryFace(f));
        }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

tf2::Field &getNodeVolume(tf2::Simulation &sim, const std::string &meshName = "")
{
    if (tf2::hasField(sim, "_NodeVolume")) return tf2::getField(sim, "_NodeVolume");

    const auto &msuf  = tf2::meshSuffix(meshName);
    const auto &domain = "Nodes"+msuf;
    const auto dim = tf2::getField(sim,"ux_N").dim;

    auto &result = tf2::newField(sim, dim, "_NodeVolume", domain);

    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);
    if (tf2::hasSMesh(sim))
    {
        const auto &m = tf2::getSMesh(sim);
        for (uint32_t it = 0; it < resultSize; it++)
        {
            const auto &N = tf2::getNodeIJKFromId(m, it);
            if (tf2::isBoundaryFace(m, N))
            {
                buffer[it] = calcFaceArea(m, N);
            }
            else
            {
                buffer[it] = tf2::calcCellVolume(m, N);
            }
        }
    }
    else
    {
        const auto &m = tf2::getUMesh(sim, meshName);
        uint32_t it = 0;
        for (auto c : m.getCellsRange())
        {
            buffer[it++] = m.getCellVolume(c);
        }
        for (auto f : m.getBFacesRange())
        {
            buffer[it++] = m.getFaceArea(f);
        }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

tf2::Field &getCellVolume(tf2::Simulation &sim, const std::string &meshName = "")
{
    if (tf2::hasField(sim, "_CellVolume")) return tf2::getField(sim, "_CellVolume");

    const auto &msuf  = tf2::meshSuffix(meshName);
    const auto &domain = "Cells"+msuf;
    const auto dim = tf2::getField(sim,"ux_N").dim;
    auto &result = tf2::newField(sim, dim, "_CellVolume", domain);

    uint32_t resultSize = tf2::getNumEntries(result);
    std::vector<double> buffer(resultSize);
    if (tf2::hasSMesh(sim))
    {
        const auto &m = tf2::getSMesh(sim);
        for (uint32_t it = 0; it < resultSize; it++)
        {
            const auto &C = tf2::getCellIJKFromId(m, it);
            buffer[it] = tf2::calcCellVolume(m, C);
        }
    }
    else
    {
        const auto &m = tf2::getUMesh(sim, meshName);
        uint32_t it = 0;
        for (auto c : m.getCellsRange())
        {
            buffer[it++] = m.getCellVolume(c);
        }
    }

    tf2::oper_setData(result, buffer.data());
    return result;
}

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
      for (uint32_t it = 0; it < resultSize; it++)
      {
        const auto &F = tf2::getFaceIJKFromId(m,it);
        buffer[it] = tf2::calcFaceStaggeredVol(m,F);
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
