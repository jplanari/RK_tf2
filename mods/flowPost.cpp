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

extern "C"
{

std::vector<double> calcFacesYPos(tf2::Simulation &sim)
{
    const auto &mesh  = tf2::getSMesh(sim);
    std::set<double, fLess> facesY;
    uint32_t NfY = TO_U32(sim.IOParamI["NY"]+1);
    const auto faces = tf2::getNumFaces(mesh);

    for (int f = 0; f < faces; ++f)
    {
        const auto &F = tf2::getFaceIJKFromId(mesh,f);
	double fy = tf2::calcFaceCentroid(mesh,F)[1];
        double nvy = tf2::calcFaceNormal(mesh,F)[1];
        if (fabs(nvy) > 0.99)
        {
            facesY.insert(fy);
            if (facesY.size() == NfY)
            {
                break;
            }
        }
    }
    TF_uAssert(facesY.size() <= NfY);

    std::vector<double> localFY;
    localFY.reserve(facesY.size());
    for (auto it : facesY)
    {
        localFY.push_back(it);
    }
    auto tmp = tf2::allGatherV(localFY);
    std::sort(tmp.begin(), tmp.end());

    std::vector<double> finalFY;
    finalFY.reserve(NfY);
    auto p = tmp.begin();
    finalFY.push_back(*p);
    p++;
    for (auto pe = tmp.end(); p != pe; p++)
    {
        if (fabs(*p-finalFY.back()) > snap)
        {
            finalFY.push_back(*p);
        }
    }

    TF_uAssert(finalFY.size() == NfY, "my size is "<<finalFY.size());
    return finalFY;
}

/*void initAvgField(tf2::Simulation &sim)
{
    auto &fab = tf2::getOrCreateField(sim, sim.IOParamI["nSims"], "avgBottom", "Faces");
    const auto &mesh = tf2::getSMesh(sim);
    std::vector<double> buf_fab(tf2::getNumFaces(mesh));
    const auto bound_faces = tf2::getNumBFaces(mesh);
    for (int bf = 0; bf<bound_faces; ++bf)
    {
        const auto &BF = tf2::getBFaceIJKFromId(mesh,bf);
        if (tf2::getBoCoId(mesh,BF) == tf2::SMesh::Bottom)
        {
            buf_fab[bf] = tf2::calcFaceArea(mesh,BF);
        }
    }
    tf2::oper_setData(fab, buf_fab.data());
    tf2::Field one(sim, sim.IOParamI["nSims"], "Faces");
    tf2::oper_setConst(one, 1.0);
    sim.IOParamD["bottomArea"] = tf2::oper_dot(one, fab)[0];
}*/


void initAvgField(tf2::Simulation &sim)
{
    auto &fab = tf2::getOrCreateField(sim, 1, "avgBottom", "Faces");
    const auto &mesh = tf2::getSMesh(sim);
    std::vector<double> buf_fab(tf2::getNumEntries(fab));
    const auto faces = tf2::getNumFaces(mesh);
    for (uint32_t it = 0; it < faces; it++)
    {
        const auto &F = tf2::getFaceIJKFromId(mesh, it);
        if (tf2::getBoCoId(mesh, F) == 3)
        {
            const auto &areas = calcCellFacesArea(mesh, F);
            buf_fab[it] = areas[1];
            // areas[1] because boco 3 is an 'y' face; for 'x' faces it would
            // have been areas[0], for 'z' faces it would have been areas[2].
        }
    }
    tf2::oper_setData(fab, buf_fab.data());
    tf2::Field one(sim, 1, "Faces");
    tf2::oper_setConst(one, 1.0);
    sim.IOParamD["bottomArea"] = tf2::oper_dot(one, fab)[0];
}

void load(tf2::Simulation &sim)
{
    const auto &name = tf2::substituteParams(sim, sim.IOParamS["load"]);
    tf2::loadFields(sim, name);
}

void averageXZ(tf2::Simulation &sim)
{
    const auto &facesYPos = calcFacesYPos(sim);
    auto NY = facesYPos.size()-1;
    std::vector<double> avguL(NY);
    std::vector<double> avgwL(NY);
    std::vector<double> avgpL(NY);
    std::vector<double> volumesL(NY);
    std::vector<double> yL(NY, -DBL_MAX);
    std::vector<double> avgdudyL(NY+1);
    std::vector<double> avgdwdyL(NY+1);
    std::vector<double> areasL(NY+1);
    std::vector<double> yfL(NY+1, -DBL_MAX);

    auto &GY = tf2::getMatrix(sim, "GradY_NF");
    auto &ux = tf2::getField(sim, "ens_avg(ux)");
    auto &uz = tf2::getField(sim, "ens_avg(uz)");
    auto &dudy = tf2::getOrCreateField(sim, 1, "dudy", "Faces");
    tf2::Field dwdy(sim, 1, "Faces");

    tf2::oper_prod(GY, ux, dudy);
    tf2::oper_prod(GY, uz, dwdy);
    auto &fab = tf2::getField(sim, "avgBottom");
    double bottomArea = sim.IOParamD["bottomArea"];
    double sdudy = fabs(tf2::oper_dot(fab, dudy)[0])/bottomArea;
    tf2::info("Calculated Re_tau is %.5e\n",sdudy);

    const auto &mesh = tf2::getSMesh(sim);
    std::vector<double> bufferu(tf2::getNumNodes(mesh));
    std::vector<double> bufferw(tf2::getNumNodes(mesh));
    std::vector<double> bufferp(tf2::getNumNodes(mesh));
    tf2::oper_getData(ux, bufferu.data());
    tf2::oper_getData(uz, bufferw.data());
    auto cells = tf2::getNumCells(mesh);
    for (uint32_t c = 0; c < cells; c++)
    {
        const auto &C = tf2::getCellIJKFromId(mesh,c);
	double cy = tf2::calcCellCentroid(mesh,C)[1];
        double vol = tf2::calcCellVolume(mesh,C);
        uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
        avguL[idp] += vol*bufferu[c];
        avgwL[idp] += vol*bufferw[c];
        volumesL[idp] += vol;
        yL[idp] = cy;
    }

    const auto &outP = sim.IOParamS["outPrefix"];
    const auto &avgu = allReduce(avguL, MPI_SUM);
    const auto &avgw = allReduce(avgwL, MPI_SUM);
    const auto &volumes = allReduce(volumesL, MPI_SUM);
    const auto &y = allReduce(yL, MPI_MAX);
    {
    std::ofstream out("results/profile_"+outP+"_cells_"+std::to_string(NY)+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ uavg wavg pavg\n";
    out<<std::scientific<<std::setprecision(8);
    if (tf2::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY; k++)
        {
            out<<y[k]
               <<" "<<y[k]*sdudy
               <<" "<<avgu[k]/volumes[k]
               <<" "<<avgw[k]/volumes[k]
               <<"\n";
        }
    }
    out.close();
    }

    std::vector<double> bufferdudy(tf2::getNumFaces(mesh));
    std::vector<double> bufferdwdy(tf2::getNumFaces(mesh));
    tf2::oper_getData(dudy, bufferdudy.data());
    tf2::oper_getData(dwdy, bufferdwdy.data());
    auto faces = tf2::getNumFaces(mesh);
    for (uint32_t f=0; f<faces; ++f)
    {
        const auto &F = tf2::getFaceIJKFromId(mesh,f);
        auto nvy = tf2::calcFaceNormal(mesh,F)[1];
        if (fabs(nvy) > 0.99)
        {
            double cy = tf2::calcFaceCentroid(mesh,F)[1];
	    const auto &areas = calcCellFacesArea(mesh,F);
            double area = fabs(areas*tf2::calcFaceNormal(mesh,F));
            uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
            avgdudyL[idp] += area*bufferdudy[f];
            avgdwdyL[idp] += area*bufferdwdy[f];
            areasL[idp] += area;
            yfL[idp] = cy;
        }
    }

    const auto &avgdudy = allReduce(avgdudyL, MPI_SUM);
    const auto &avgdwdy = allReduce(avgdwdyL, MPI_SUM);
    const auto &areas = allReduce(areasL, MPI_SUM);
    const auto &yf = allReduce(yfL, MPI_MAX);
    {
    std::ofstream out("results/profile_"+outP+"_faces_"+std::to_string(NY)+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ duavg/dy dwavg/dy\n";
    out<<std::scientific<<std::setprecision(8);
    if (tf2::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY+1; k++)
        {
            out<<yf[k]
               <<" "<<yf[k]*sdudy
               <<" "<<avgdudy[k]/areas[k]
               <<" "<<avgdwdy[k]/areas[k]
               <<"\n";
        }
    }
    out.close();
    }
}

void ReynoldsStresses(tf2::Simulation &sim)
{
    auto &avg_u = tf2::getField(sim, "ens_avg(ux)");
    auto &avg_v = tf2::getField(sim, "ens_avg(uy)");
    auto &avg_w = tf2::getField(sim, "ens_avg(uz)");
    auto &avg_uu = tf2::getField(sim, "ens_avg(ux*ux)");
    auto &avg_uv = tf2::getField(sim, "ens_avg(ux*uy)");
    auto &avg_uw = tf2::getField(sim, "ens_avg(ux*uz)");
    auto &avg_vv = tf2::getField(sim, "ens_avg(uy*uy)");
    auto &avg_vw = tf2::getField(sim, "ens_avg(uy*uz)");
    auto &avg_ww = tf2::getField(sim, "ens_avg(uz*uz)");

    tf2::Field R_uu(sim, 1, "Nodes");
    tf2::Field R_uv(sim, 1, "Nodes");
    tf2::Field R_uw(sim, 1, "Nodes");
    tf2::Field R_vv(sim, 1, "Nodes");
    tf2::Field R_vw(sim, 1, "Nodes");
    tf2::Field R_ww(sim, 1, "Nodes");

    tf2::oper_fmadd(avg_u, avg_u, R_uu, 1.0, 0.0);
    tf2::oper_axpy(avg_uu, R_uu, 1.0, -1.0);

    tf2::oper_fmadd(avg_u, avg_v, R_uv, 1.0, 0.0);
    tf2::oper_axpy(avg_uv, R_uv, 1.0, -1.0);

    tf2::oper_fmadd(avg_u, avg_w, R_uw, 1.0, 0.0);
    tf2::oper_axpy(avg_uw, R_uw, 1.0, -1.0);

    tf2::oper_fmadd(avg_v, avg_v, R_vv, 1.0, 0.0);
    tf2::oper_axpy(avg_vv, R_vv, 1.0, -1.0);

    tf2::oper_fmadd(avg_v, avg_w, R_vw, 1.0, 0.0);
    tf2::oper_axpy(avg_vw, R_vw, 1.0, -1.0);

    tf2::oper_fmadd(avg_w, avg_w, R_ww, 1.0, 0.0);
    tf2::oper_axpy(avg_ww, R_ww, 1.0, -1.0);

    const auto &facesYPos = calcFacesYPos(sim);
    auto NY = facesYPos.size()-1;
    std::vector<double> avgLR_uu(NY);
    std::vector<double> avgLR_uv(NY);
    std::vector<double> avgLR_uw(NY);
    std::vector<double> avgLR_vv(NY);
    std::vector<double> avgLR_vw(NY);
    std::vector<double> avgLR_ww(NY);
    std::vector<double> volumesL(NY);
    std::vector<double> yL(NY, -DBL_MAX);
    const auto &mesh = tf2::getSMesh(sim);
    std::vector<double> buffer_uu(tf2::getNumNodes(mesh));
    std::vector<double> buffer_uv(tf2::getNumNodes(mesh));
    std::vector<double> buffer_uw(tf2::getNumNodes(mesh));
    std::vector<double> buffer_vv(tf2::getNumNodes(mesh));
    std::vector<double> buffer_vw(tf2::getNumNodes(mesh));
    std::vector<double> buffer_ww(tf2::getNumNodes(mesh));
    tf2::oper_getData(R_uu, buffer_uu.data());
    tf2::oper_getData(R_uv, buffer_uv.data());
    tf2::oper_getData(R_uw, buffer_uw.data());
    tf2::oper_getData(R_vv, buffer_vv.data());
    tf2::oper_getData(R_vw, buffer_vw.data());
    tf2::oper_getData(R_ww, buffer_ww.data());
    auto const cells = tf2::getNumCells(mesh);
    for (int c = 0; c<cells; ++c)
    {
        const auto &C = tf2::getCellIJKFromId(mesh,c);
        double cy = tf2::calcCellCentroid(mesh,C)[1];
        double vol = tf2::calcCellVolume(mesh,C);
        uint64_t idp = bsearchVectorE(facesYPos, cy, snap);
        avgLR_uu[idp] += vol*buffer_uu[c];
        avgLR_uv[idp] += vol*buffer_uv[c];
        avgLR_uw[idp] += vol*buffer_uw[c];
        avgLR_vv[idp] += vol*buffer_vv[c];
        avgLR_vw[idp] += vol*buffer_vw[c];
        avgLR_ww[idp] += vol*buffer_ww[c];
        volumesL[idp] += vol;
        yL[idp] = cy;
    }

    auto &GY = tf2::getMatrix(sim, "GradY_NF");
    tf2::Field dudy(sim, 1, "Faces");
    tf2::oper_prod(GY, avg_u, dudy);
    auto &fab = tf2::getField(sim, "avgBottom");
    double bottomArea = sim.IOParamD["bottomArea"];
    double sdudy = fabs(tf2::oper_dot(fab, dudy)[0])/bottomArea;

    const auto &outP = sim.IOParamS["outPrefix"];
    const auto &avgR_uu = allReduce(avgLR_uu, MPI_SUM);
    const auto &avgR_uv = allReduce(avgLR_uv, MPI_SUM);
    const auto &avgR_uw = allReduce(avgLR_uw, MPI_SUM);
    const auto &avgR_vv = allReduce(avgLR_vv, MPI_SUM);
    const auto &avgR_vw = allReduce(avgLR_vw, MPI_SUM);
    const auto &avgR_ww = allReduce(avgLR_ww, MPI_SUM);
    const auto &volumes = allReduce(volumesL, MPI_SUM);
    const auto &y = allReduce(yL, MPI_MAX);
    {
    std::ofstream out("results/reystress_"+outP+"_"+std::to_string(NY)+".dat");
    TF_uAssert(out.is_open(), "could not open filename for output");
    out<<"#y y+ R_uu R_vv R_ww R_uv R_uw R_vw\n";
    out<<std::scientific<<std::setprecision(8);
    if (tf2::mpiRank() == 0)
    {
        for (uint32_t k = 0; k < NY; k++)
        {
            out<<y[k]
               <<" "<<y[k]*sdudy
               <<" "<<avgR_uu[k]/volumes[k]
               <<" "<<avgR_vv[k]/volumes[k]
               <<" "<<avgR_ww[k]/volumes[k]
               <<" "<<avgR_uv[k]/volumes[k]
               <<" "<<avgR_uw[k]/volumes[k]
               <<" "<<avgR_vw[k]/volumes[k]
               <<"\n";
        }
    }
    out.close();
    }
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
    auto &uxux = tf2::getField(sim, "avg(ux_N*ux_N)");
    auto &uxuy = tf2::getField(sim, "avg(ux_N*uy_N)");
    auto &uxuz = tf2::getField(sim, "avg(ux_N*uz_N)");
    auto &uyuy = tf2::getField(sim, "avg(uy_N*uy_N)");
    auto &uyuz = tf2::getField(sim, "avg(uy_N*uz_N)");
    auto &uzuz = tf2::getField(sim, "avg(uz_N*uz_N)");



    auto &ensAvg_ux = tf2::getOrCreateField(sim, 1, "ens_avg(ux)", "Nodes");
    auto &ensAvg_uy = tf2::getOrCreateField(sim, 1, "ens_avg(uy)", "Nodes");
    auto &ensAvg_uz = tf2::getOrCreateField(sim, 1, "ens_avg(uz)", "Nodes");
    auto &ensAvg_uxux = tf2::getOrCreateField(sim, 1, "ens_avg(ux*ux)", "Nodes");
    auto &ensAvg_uxuy = tf2::getOrCreateField(sim, 1, "ens_avg(ux*uy)", "Nodes");
    auto &ensAvg_uxuz = tf2::getOrCreateField(sim, 1, "ens_avg(ux*uz)", "Nodes");
    auto &ensAvg_uyuy = tf2::getOrCreateField(sim, 1, "ens_avg(uy*uy)", "Nodes");
    auto &ensAvg_uyuz = tf2::getOrCreateField(sim, 1, "ens_avg(uy*uz)", "Nodes");
    auto &ensAvg_uzuz = tf2::getOrCreateField(sim, 1, "ens_avg(uz*uz)", "Nodes");

    apply_ensemble_average(ux, ensAvg_ux);
    apply_ensemble_average(uy, ensAvg_uy);
    apply_ensemble_average(uz, ensAvg_uz);
    apply_ensemble_average(uxux, ensAvg_uxux);
    apply_ensemble_average(uxuy, ensAvg_uxuy);
    apply_ensemble_average(uxuz, ensAvg_uxuz);
    apply_ensemble_average(uyuy, ensAvg_uyuy);
    apply_ensemble_average(uyuz, ensAvg_uyuz);
    apply_ensemble_average(uzuz, ensAvg_uzuz);
}

}
