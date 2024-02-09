#include <algorithm>
#include <cfloat>
#include <iomanip>
#include "tf2/DataIO.h"
#include "tf2/IO.h"
#include "tf2/Simulation.h"
#include "tf2/Opers.h"
#include "tf2/ll_Opers.h"

static const double snap = 1e-8;

struct fLess
{
    constexpr bool operator()(const double &lhs, const double &rhs) const
    {
        return (lhs < rhs && fabs(lhs-rhs) > snap);
    }
};

uint64_t bsearchVectorE(const std::vector<double> &vec, double target, double eps)
{
    // Given a vector 'vec' such that vec[i]>vec[i-1], find the largest i such
    // that vec[i]<=target using a binary search.
    uint64_t li = 0, ui = vec.size()-1;
    while ((ui-li)>1)
    {
        uint64_t mi = (li+ui)/2;
        if (fabs(vec[mi]-target) < eps)
        {
            return mi;
        }
        else if (vec[mi]>target)
        {
            ui = mi;
        }
        else
        {
            li = mi;
        }
    }
    return vec[ui] <= target? ui : li;
}

std::vector<double> allReduce(const std::vector<double> &in, MPI_Op op, MPI_Comm w = MPI_COMM_WORLD)
{
    std::vector<double> result(in.size());
    int32_t sz = TO_S32(in.size());
    tf2::checkr(MPI_Allreduce(in.data(), result.data(), sz, MPI_DOUBLE, op, w),
                "allreduce");
    return result;
}

TF_Func void SetUp_Matrices(tf2::Simulation &sim)
{
    // This function exists basically to initialize all the fields required by
    // the official FSM, with the correct number of components, i.e., one
    // component per simulation.
    const auto &meshName = tf2::getMeshName(sim);

    const auto smesh = tf2::hasSMesh(sim);
   
    const auto &pMod = smesh? "smeshPatterns.so" : "patterns.so";
    const auto &kMod = smesh? "smeshInterpolators.so" : "interpolators.so";
    
    if (not tf2::hasMatrix(sim, "Id_NC"))
    {
        tf2::newMatrix(sim, "Id_NC", "pId_NC", pMod, "kId_NC", kMod, meshName);
    }
    if (not tf2::hasMatrix(sim, "Id_CN"))
    {
        tf2::newMatrix(sim, "Id_CN", "pId_CN", pMod, "kId_CN", kMod, meshName);
    }
}

extern "C"
{

std::vector<double> calcFacesZPos(tf2::Simulation &sim)
{
    const auto &mesh  = tf2::getSMesh(sim);
    std::set<double, fLess> facesZ;
    uint32_t NfZ = TO_U32(sim.IOParamI["NZ"]+1);
    const auto faces = tf2::getNumFaces(mesh);

    for (int f = 0; f < faces; ++f)
    {
        const auto &F = tf2::getFaceIJKFromId(mesh,f);
	      double fz = tf2::calcFaceCentroid(mesh,F)[2];
        double nvz = tf2::calcFaceNormal(mesh,F)[2];
        if (fabs(nvz) > 0.99)
        {
            facesZ.insert(fz);
            if (facesZ.size() == NfZ)
            {
                break;
            }
        }
    }
    TF_uAssert(facesZ.size() <= NfZ);

    std::vector<double> localFZ;
    localFZ.reserve(facesZ.size());
    for (auto it : facesZ)
    {
        localFZ.push_back(it);
    }
    auto tmp = tf2::allGatherV(localFZ);
    std::sort(tmp.begin(), tmp.end());

    std::vector<double> finalFZ;
    finalFZ.reserve(NfZ);
    auto p = tmp.begin();
    finalFZ.push_back(*p);
    p++;
    for (auto pe = tmp.end(); p != pe; p++)
    {
        if (fabs(*p-finalFZ.back()) > snap)
        {
            finalFZ.push_back(*p);
        }
    }

    TF_uAssert(finalFZ.size() == NfZ, "my size is "<<finalFZ.size());
    return finalFZ;
}

void initAvgField(tf2::Simulation &sim)
{
    auto &fab = tf2::getOrCreateField(sim, 1, "avgBack", "Faces");
    const auto &mesh = tf2::getSMesh(sim);
    std::vector<double> buf_fab(tf2::getNumEntries(fab));
    const auto faces = tf2::getNumFaces(mesh);
    for (uint32_t it = 0; it < faces; it++)
    {
        const auto &F = tf2::getFaceIJKFromId(mesh, it);
        if (tf2::getBoCoId(mesh, F) == 4)
        {
            const auto &areas = calcCellFacesArea(mesh, F);
            buf_fab[it] = areas[2];
            // areas[1] because boco 3 is an 'y' face; for 'x' faces it would
            // have been areas[0], for 'z' faces it would have been areas[2].
        }
    }
    tf2::oper_setData(fab, buf_fab.data());
    tf2::Field one(sim, 1, "Faces");
    tf2::oper_setConst(one, 1.0);
    sim.IOParamD["backArea"] = tf2::oper_dot(one, fab)[0];
}



void load(tf2::Simulation &sim)
{
    const auto &name = tf2::substituteParams(sim, sim.IOParamS["load"]);
    tf2::loadFields(sim, name);
}

void apply_ensemble_average(tf2::Field &src, tf2::Field &dest)
{
    TF_uAssert(dest.dim == 1);
    TF_uAssert(&tf2::getDomain(src) == &tf2::getDomain(dest));

    std::vector<double> srcBuffer(tf2::getNumEntries(src));
    std::vector<double> destBuffer(tf2::getNumEntries(dest));

    tf2::oper_getData(src, srcBuffer.data());

    auto srcElem = tf2::getDomain(src).localSize;
    for (uint32_t k = 0; k < srcElem; k++)
    {
        double avg = 0.0;
        for (uint32_t d = 0; d < src.dim; d++) avg += srcBuffer[k*src.dim+d];
        destBuffer[k] = avg/TO_F64(src.dim);
    }

    tf2::oper_setData(dest, destBuffer.data());
}

void ensemble_avg(tf2::Simulation &sim)
{
    auto &ux = tf2::getField(sim, "avg(ux_N)");
    auto &uy = tf2::getField(sim, "avg(uy_N)");
    auto &uz = tf2::getField(sim, "avg(uz_N)");
    auto &T = tf2::getField(sim, "avg(T_N)");

    auto &ensAvg_ux = tf2::getOrCreateField(sim, 1, "ens_avg(ux)", "Nodes");
    auto &ensAvg_uy = tf2::getOrCreateField(sim, 1, "ens_avg(uy)", "Nodes");
    auto &ensAvg_uz = tf2::getOrCreateField(sim, 1, "ens_avg(uz)", "Nodes");
    auto &ensAvg_T = tf2::getOrCreateField(sim, 1, "ens_avg(T)", "Nodes");

    apply_ensemble_average(ux, ensAvg_ux);
    apply_ensemble_average(uy, ensAvg_uy);
    apply_ensemble_average(uz, ensAvg_uz);
    apply_ensemble_average(T, ensAvg_T);
}

void computeNusselt(tf2::Simulation &sim)
{
  auto &T = tf2::getField(sim, "ens_avg(T)");
  auto &GX = tf2::getMatrix(sim, "GradX_CC");
  auto &ID_NC = tf2::getMatrix(sim, "Id_NC");

  auto NX = sim.IOParamI["NX"];
  auto &Tc = TF_getTmpField(sim, 1, "Cells");
  tf2::oper_prod(ID_NC, T, Tc);

  auto &dTdx = TF_getTmpField(sim, Tc);
  tf2::oper_prod(GX, Tc, dTdx);
  
  const auto &m = tf2::getSMesh(sim);
  std::vector<double> dTdx_buf(tf2::getNumEntries(dTdx));
  std::vector<double> Nux_buf(tf2::getNumEntries(dTdx));
  std::vector<double> left_buf(tf2::getNumEntries(dTdx));
  std::vector<double> right_buf(tf2::getNumEntries(dTdx));
  tf2::oper_getData(dTdx, dTdx_buf.data());

  auto srcElem = tf2::getDomain(dTdx).localSize;
  double dy,dz;
  for (uint32_t k=0; k < srcElem; k++)
  {
    auto C = tf2::getCellIJKFromId(m,k);

    auto Ftop = tf2::getFaceIJKFromId(m,tf2::getFaceIndex(m,C.i,C.j,C.k,3));
    auto Fbot = tf2::getFaceIJKFromId(m,tf2::getFaceIndex(m,C.i,C.j,C.k,2));
    
    dz = tf2::calcFaceCentroid(m,Ftop)[2]-tf2::calcFaceCentroid(m,Fbot)[2];

    auto Ffro = tf2::getFaceIJKFromId(m,tf2::getFaceIndex(m,C.i,C.j,C.k,5));
    auto Fbac = tf2::getFaceIJKFromId(m,tf2::getFaceIndex(m,C.i,C.j,C.k,4));

    dy = tf2::calcFaceCentroid(m,Ffro)[1] - tf2::calcFaceCentroid(m,Fbac)[1];

    Nux_buf[k] = dTdx_buf[k]*dy*dz;

    if (C.i==0) left_buf[k] = 1.0;
    else left_buf[k] = 0.0;

    if (C.i==NX) right_buf[k] = 1.0;
    else right_buf[k] = 0.0;
  }
  
  auto &Nu = TF_getTmpField(sim, dTdx);
  tf2::oper_setData(Nu, Nux_buf.data());
  
  auto &left = TF_getTmpField(sim, Nu);
  tf2::oper_setData(left,left_buf.data());
  tf2::info("Left wall Nusselt = %e\n",tf2::oper_dot(Nu,left)[0]/4);

  auto &right = TF_getTmpField(sim, Nu);
  tf2::oper_setData(right,right_buf.data());
  tf2::info("Right wall Nusselt = %e\n",tf2::oper_dot(Nu,right)[0]/4);
}
}
