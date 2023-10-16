void printMaxVals(std::string point, tf2::Field &ux, tf2::Field &uy, tf2::Field &uz)
{
    tf2::info("max values at this point: %s\n",point.c_str());
    tf2::info("\tux=%e\n",tf2::max(tf2::oper_max(ux)));
    tf2::info("\tuy=%e\n",tf2::max(tf2::oper_max(uy)));
    tf2::info("\tuz=%e\n",tf2::max(tf2::oper_max(uz)));
}

static void calcStaggeredVolumes(const tf2::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getFacesRange())
    {
        const auto &nbs = mesh.getFaceNbCells(it);
        const auto &fnv = mesh.getFaceNormal(it);
        tf2::V3 d;
        if (nbs.second != mt::CellIdT::None)
        {
            d = mesh.getCellCentroid(nbs.second)-mesh.getCellCentroid(nbs.first);
        }
        else
        {
            d = mesh.getFaceCentroid(it)-mesh.getCellCentroid(nbs.first);
        }
        *buffer++ = (mesh.getFaceArea(it)*fabs(fnv*d));
    }
}

static void calcVolumes(const tf2::SMesh &mesh, double *buffer)
{
    const auto nc = tf2::getNumCells(mesh);
    for (uint32_t it = 0; it < nc; it++)
    {
        const auto C = tf2::getCellIJKFromId(mesh, it);
        *buffer++ = tf2::calcCellVolume(mesh, C);
    }
}

static void calcVolumes(const tf2::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getCellsRange())
    {
        *buffer++ = mesh.getCellVolume(it);
    }
}

static void calcInvVolumes(const tf2::UMesh &mesh, double *buffer)
{
    for (auto it : mesh.getCellsRange())
    {
        *buffer++ = 1.0/mesh.getCellVolume(it);
    }
}

static void genVolumesField(tf2::Simulation &sim)
{
    if (not tf2::hasField(sim, "Omega_C"))
    {
        std::vector<double> buffer(tf2::getDomainSize(sim, "Cells"));
        if (hasSMesh(sim))
        {
            const auto &m = tf2::getSMesh(sim);
            calcVolumes(m, buffer.data());
        }
        else
        {
            const auto &m = tf2::getUMesh(sim);
            calcVolumes(m, buffer.data());
        }
        auto &f = tf2::newField(sim, 1, "Omega_C", "Cells");
        tf2::oper_setData(f, buffer.data());
    }
}

TF_Func void SetUp_Momentum_Collocated(tf2::Simulation &sim)
{
    TF_uAssert((sim.meshes.size() == 1) xor tf2::hasSMesh(sim), "not yet implemented for multiple meshes");

    const auto &meshName = tf2::getMeshName(sim);
    const auto &msuf     = tf2::meshSuffix(meshName);
    const auto &cells    = "Cells"+msuf;
    const auto &faces    = "Faces"+msuf;
    const auto &nodes    = "Nodes"+msuf;

    const auto smesh = tf2::hasSMesh(sim);
    const auto &op   = tf2::getOfficialPath();
    const auto &pMod = op+(smesh? "/smeshPatterns.so" : "/patterns.so");
    const auto &kMod = op+(smesh? "/smeshInterpolators.so" : "/interpolators.so");
    if (not tf2::hasMatrix(sim, "Id_NC"))
    {
        tf2::newMatrix(sim, "Id_NC", "pId_NC", pMod, "kId_NC", kMod, meshName);
    }
    if (not tf2::hasMatrix(sim, "Id_CN"))
    {
        tf2::newMatrix(sim, "Id_CN", "pId_CN", pMod, "kId_CN", kMod, meshName);
    }

    genVolumesField(sim);

    tf2::getOrCreateField(sim, 1, "upredx_C",  cells);
    tf2::getOrCreateField(sim, 1, "upredy_C",  cells);
    tf2::getOrCreateField(sim, 1, "upredz_C",  cells);
    tf2::getOrCreateField(sim, 1, "momSrcx_C", cells);
    tf2::getOrCreateField(sim, 1, "momSrcy_C", cells);
    tf2::getOrCreateField(sim, 1, "momSrcz_C", cells);
    tf2::getOrCreateField(sim, 1, "P_C",       cells);

    tf2::getOrCreateField(sim, 1, "ux_F", faces);
    tf2::getOrCreateField(sim, 1, "uy_F", faces);
    tf2::getOrCreateField(sim, 1, "uz_F", faces);

 }


