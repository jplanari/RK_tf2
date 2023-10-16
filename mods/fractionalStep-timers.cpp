#include "tf2/ll_Opers.h"
#include "tf2/Opers.h"
#include "tf2/Simulation.h"
#include "butcherTableaux.h"

#include <vector>
#include <cstdarg>

// +---------------------------------------------------------------------------+
// |                              General remarks                              |
// +---------------------------------------------------------------------------+

// This module introduces some functions to help solving the Navier-Stokes
// equations. For the moment only incompressible fluids have been considered.

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

TF_Func void SetUp_Momentum_Staggered(tf2::Simulation &sim)
{
    
    SetUp_Momentum_Collocated(sim);

    const bool hasInterpCF = tf2::hasMatrix(sim, "Interp_CF");
    const bool hasInterpFC = tf2::hasMatrix(sim, "Interp_FC");

    if (hasInterpCF && hasInterpFC) return;

    const auto &meshName = tf2::getMeshName(sim);
    const auto &msuf     = tf2::meshSuffix(meshName);
    const auto &cells    = "Cells"+msuf;
    const auto &faces    = "Faces"+msuf;

    const auto smesh = tf2::hasSMesh(sim);
    const auto &op   = tf2::getOfficialPath();
    const auto &pMod = op+(smesh? "/smeshPatterns.so" : "/patterns.so");
    const auto &kMod = op+(smesh? "/smeshInterpolators.so" : "/interpolators.so");
    if (not hasInterpCF && not hasInterpFC)
    {
        tf2::newMatrix(sim, "Interp_CF", "p_CF_FirstNb", pMod, "k_Interp_CF_VW_ZB_FirstNb", kMod, meshName);
        tf2::newMatrix(sim, "Interp_FC", "p_FC_FirstNb", pMod, "k_Interp_FC_VW_ZB_FirstNb", kMod, meshName);
    }
    else if (not hasInterpFC)
    {
        const auto &mesh = tf2::getUMesh(sim);
        tf2::Field InvOmega_C(sim, 1, cells);
        std::vector<double> bc(tf2::getNumEntries(InvOmega_C));
        calcInvVolumes(mesh, bc.data());
        tf2::oper_setData(InvOmega_C, bc.data());

        tf2::Field Omega_F(sim, 1, faces);
        std::vector<double> bf(tf2::getNumEntries(Omega_F));
        calcStaggeredVolumes(mesh, bf.data());
        tf2::oper_setData(Omega_F, bf.data());

        const auto &ICF = tf2::getMatrix(sim, "Interp_CF");
        auto mt = tf2::ll_transpose(ICF.csr, tf2::getDomainSize(sim, cells));
        auto p1 = tf2::ll_prodDiag(mt, Omega_F.array);
        auto p2 = tf2::ll_prodDiag(InvOmega_C.array, p1);
        tf2::newMatrix(sim, "Interp_FC", faces, cells, p2);
    }
    else
    {
        const auto &mesh = tf2::getUMesh(sim);
        tf2::Field Omega_F(sim, 1, faces);
        std::vector<double> bf(tf2::getNumEntries(Omega_F));
        calcStaggeredVolumes(mesh, bf.data());
        tf2::oper_setData(Omega_F, bf.data());
        tf2::Field InvOmega_F(sim, 1, faces);
        tf2::oper_setConst(InvOmega_F, 1.0);
        tf2::oper_div(InvOmega_F, Omega_F);

        tf2::Field InvOmega_C(sim, 1, cells);
        std::vector<double> bc(tf2::getNumEntries(InvOmega_C));
        calcInvVolumes(mesh, bc.data());
        tf2::oper_setData(InvOmega_C, bc.data());
        tf2::Field Omega_C(sim, 1, cells);
        tf2::oper_setConst(Omega_C, 1.0);
        tf2::oper_div(Omega_C, InvOmega_C);

        const auto &IFC = tf2::getMatrix(sim, "Interp_FC");
        auto mt = tf2::ll_transpose(IFC.csr, tf2::getDomainSize(sim, faces));
        auto p1 = tf2::ll_prodDiag(mt, Omega_C.array);
        auto p2 = tf2::ll_prodDiag(InvOmega_F.array, p1);
        tf2::newMatrix(sim, "Interp_CF", faces, cells, p2);
    }
}

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

tf2::Field &predictor(tf2::Field &un, tf2::Field &conv, tf2::Field &diff, tf2::Field &ms, double coef, double dt, tf2::Simulation &sim)
{
    auto &up = TF_getTmpField(sim,un);

    tf2::oper_axpy(un, up, 1.0, 0.0);
    tf2::oper_axpy(diff, up, dt*coef, 1.0);
    tf2::oper_axpy(conv, up, -dt*coef, 1.0);
    tf2::oper_axpy(ms, up, dt*coef, 1.0);

    return up;
}

tf2::Field &predictorVector(tf2::Field &un, std::vector<tf2::Field*> conv, std::vector<tf2::Field*> diff, tf2::Field &ms, butcherTableau coefs, double dt, int i, tf2::Simulation &sim)
{
  auto &up = TF_getTmpField(sim,un);

  if(i!=0)
  {
    up = predictor(un,*(conv.at(0)),*(diff.at(0)),ms,coefs.A.at(id(i,1)),dt,sim);
    for (int j=2; j<i; ++j)
      up = predictor(up,*(conv.at(j-1)),*(diff.at(j-1)),ms,coefs.A.at(id(i,j)),dt,sim);
    return up;
  }
  else
  {
    up = predictor(un,*(conv.at(0)),*(diff.at(0)),ms,coefs.b.at(0),dt,sim);
    for (int j=1; j<coefs.b.size(); ++j)
      up = predictor(up,*(conv.at(j)),*(diff.at(j)),ms,coefs.b.at(j),dt,sim);
    return up;
  }
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

void printMaxVals(std::string point, tf2::Field &ux, tf2::Field &uy, tf2::Field &uz)
{
    tf2::info("max values at this point: %s\n",point.c_str());
    tf2::info("\tux=%e\n",tf2::max(tf2::oper_max(ux)));
    tf2::info("\tuy=%e\n",tf2::max(tf2::oper_max(uy)));
    tf2::info("\tuz=%e\n",tf2::max(tf2::oper_max(uz)));
}

TF_Func bool EulerIteration(tf2::Simulation &sim)
{
    auto &M     = tf2::getMatrix(sim, "Id_NC");
    auto &GX    = tf2::getMatrix(sim, "GradX_CF");
    auto &GY    = tf2::getMatrix(sim, "GradY_CF");
    auto &GZ    = tf2::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tf2::getMatrix(sim, "Interp_CF");
    auto &IFC   = tf2::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
    
    auto &ux     = tf2::getField(sim, "ux_N");
    auto &uy     = tf2::getField(sim, "uy_N");
    auto &uz     = tf2::getField(sim, "uz_N");
    
    auto &ufx    = tf2::getField(sim, "ux_F");
    auto &ufy    = tf2::getField(sim, "uy_F");
    auto &ufz    = tf2::getField(sim, "uz_F");

    auto &upx    = tf2::getField(sim, "upredx_C");
    auto &upy    = tf2::getField(sim, "upredy_C");
    auto &upz    = tf2::getField(sim, "upredz_C");
    auto &msx    = tf2::getField(sim, "momSrcx_C");
    auto &msy    = tf2::getField(sim, "momSrcy_C");
    auto &msz    = tf2::getField(sim, "momSrcz_C");
    
    auto &p      = tf2::getField(sim, "P_C");
    static auto psolver = tf2::getSolver(sim, "Pressure_Solver");
   
    // Temporary variables.
    auto &pSource = TF_getTmpField(sim, p);
    
    auto &uxn = TF_getTmpField(sim,upx);
    auto &uyn = TF_getTmpField(sim,upx);
    auto &uzn = TF_getTmpField(sim,upx);

    tf2::oper_prod(M, ux, upx);
    tf2::oper_prod(M, uy, upy);
    tf2::oper_prod(M, uz, upz);

    tf2::oper_copy(upx,uxn);
    tf2::oper_copy(upy,uyn);
    tf2::oper_copy(upz,uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    //printMaxVals("end of first stage",ux,uy,uz);

    //FINAL STAGE
    
    auto &diffx = diffusive(ux,0,sim);
    auto &diffy = diffusive(uy,1,sim);
    auto &diffz = diffusive(uz,2,sim);
    auto &convx = convective(ufx,0,sim);
    auto &convy = convective(ufy,1,sim);
    auto &convz = convective(ufz,2,sim);

    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    
    upx = predictor(uxn,convx,diffx,msx,dt,1.0,sim);
    upy = predictor(uyn,convy,diffy,msy,dt,1.0,sim);
    upz = predictor(uzn,convz,diffz,msz,dt,1.0,sim);
   
    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tf2::oper_solve(psolver, pSource, p);
 
    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
   // Map the velocity to the nodes.
    tf2::oper_prod(ID_CN, upx, ux);
    tf2::oper_prod(ID_CN, upy, uy);
    tf2::oper_prod(ID_CN, upz, uz);

    //printMaxVals("end of time-step",ux,uy,uz);

    return tf2::Iter_Continue;
}

TF_Func bool RK2Iteration(tf2::Simulation &sim)
{
    auto &M     = tf2::getMatrix(sim, "Id_NC");
    auto &GX    = tf2::getMatrix(sim, "GradX_CF");
    auto &GY    = tf2::getMatrix(sim, "GradY_CF");
    auto &GZ    = tf2::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tf2::getMatrix(sim, "Interp_CF");
    auto &IFC   = tf2::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
    
    auto &ux     = tf2::getField(sim, "ux_N");
    auto &uy     = tf2::getField(sim, "uy_N");
    auto &uz     = tf2::getField(sim, "uz_N");
    
    auto &ufx    = tf2::getField(sim, "ux_F");
    auto &ufy    = tf2::getField(sim, "uy_F");
    auto &ufz    = tf2::getField(sim, "uz_F");

    auto &upx    = tf2::getField(sim, "upredx_C");
    auto &upy    = tf2::getField(sim, "upredy_C");
    auto &upz    = tf2::getField(sim, "upredz_C");

    auto &msx    = tf2::getField(sim, "momSrcx_C");
    auto &msy    = tf2::getField(sim, "momSrcy_C");
    auto &msz    = tf2::getField(sim, "momSrcz_C");
    
    auto &p      = tf2::getField(sim, "P_C");
    static auto psolver = tf2::getSolver(sim, "Pressure_Solver");

   
    // Temporary variables.
    auto &pSource = TF_getTmpField(sim, p);
    
    auto &uxn = TF_getTmpField(sim,upx);
    auto &uyn = TF_getTmpField(sim,upx);
    auto &uzn = TF_getTmpField(sim,upx);

    tf2::oper_prod(M, ux, upx);
    tf2::oper_prod(M, uy, upy);
    tf2::oper_prod(M, uz, upz);

    tf2::oper_copy(upx,uxn);
    tf2::oper_copy(upy,uyn);
    tf2::oper_copy(upz,uzn);

    double dt = sim.IOParamD["_TimeStep"];

    //SUBSTAGE 1
    
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    //SUBSTAGE 2
    auto &diffx  = diffusive(ux,0,sim);
    auto &diffy  = diffusive(uy,1,sim);
    auto &diffz  = diffusive(uz,2,sim);
    auto &convx  = convective(ufx,0,sim);
    auto &convy  = convective(ufy,1,sim);
    auto &convz  = convective(ufz,2,sim);

    //predictor 2 stage

    upx = predictor(uxn,convx,diffx,msx,1.0,dt,sim);
    upy = predictor(uyn,convy,diffy,msy,1.0,dt,sim);
    upz = predictor(uzn,convz,diffz,msz,1.0,dt,sim);
   
    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tf2::oper_solve(psolver, pSource, p);
 
    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
   // Map the velocity to the nodes.
    tf2::oper_prod(ID_CN, upx, ux);
    tf2::oper_prod(ID_CN, upy, uy);
    tf2::oper_prod(ID_CN, upz, uz);

    //FINAL STAGE
 
    auto &diffx2  = diffusive(ux,0,sim);
    auto &diffy2  = diffusive(uy,1,sim);
    auto &diffz2  = diffusive(uz,2,sim);
    auto &convx2  = convective(ufx,0,sim);
    auto &convy2  = convective(ufy,1,sim);
    auto &convz2  = convective(ufz,2,sim);
    
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    
    upx = predictor(uxn,convx,diffx,msx,0.5,dt,sim);
    upy = predictor(uyn,convy,diffy,msy,0.5,dt,sim);
    upz = predictor(uzn,convz,diffz,msz,0.5,dt,sim);

    upx = predictor(upx,convx2,diffx2,msx,0.5,dt,sim);
    upy = predictor(upy,convy2,diffy2,msy,0.5,dt,sim);
    upz = predictor(upz,convz2,diffz2,msz,0.5,dt,sim);

    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tf2::oper_solve(psolver, pSource, p);
 
    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
   // Map the velocity to the nodes.
    tf2::oper_prod(ID_CN, upx, ux);
    tf2::oper_prod(ID_CN, upy, uy);
    tf2::oper_prod(ID_CN, upz, uz);

    //printMaxVals("end of time-step",ux,uy,uz);

    return tf2::Iter_Continue;
}

TF_Func bool RKiteration(tf2::Simulation &sim)
{
    auto &M     = tf2::getMatrix(sim, "Id_NC");
    auto &GX    = tf2::getMatrix(sim, "GradX_CF");
    auto &GY    = tf2::getMatrix(sim, "GradY_CF");
    auto &GZ    = tf2::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tf2::getMatrix(sim, "Interp_CF");
    auto &IFC   = tf2::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
    
    auto &ux     = tf2::getField(sim, "ux_N");
    auto &uy     = tf2::getField(sim, "uy_N");
    auto &uz     = tf2::getField(sim, "uz_N");
    
    auto &ufx    = tf2::getField(sim, "ux_F");
    auto &ufy    = tf2::getField(sim, "uy_F");
    auto &ufz    = tf2::getField(sim, "uz_F");

    auto &upx    = tf2::getField(sim, "upredx_C");
    auto &upy    = tf2::getField(sim, "upredy_C");
    auto &upz    = tf2::getField(sim, "upredz_C");

    auto &msx    = tf2::getField(sim, "momSrcx_C");
    auto &msy    = tf2::getField(sim, "momSrcy_C");
    auto &msz    = tf2::getField(sim, "momSrcz_C");
    
    auto &p      = tf2::getField(sim, "P_C");
    static auto psolver = tf2::getSolver(sim, "Pressure_Solver");
    
    butcherTableau coefs = intScheme("heunRK2");

    tf2::info("Number of stages per iteration: %d\n",coefs.b.size());

    double dt = sim.IOParamD["_TimeStep"];
    
    std::vector<tf2::Field*> convx(coefs.b.size());
    std::vector<tf2::Field*> diffx(coefs.b.size());
    std::vector<tf2::Field*> convy(coefs.b.size());
    std::vector<tf2::Field*> diffy(coefs.b.size());
    std::vector<tf2::Field*> convz(coefs.b.size());
    std::vector<tf2::Field*> diffz(coefs.b.size());
    
    // Temporary variables.
    auto &pSource = TF_getTmpField(sim, p);
    
    auto &uxn = TF_getTmpField(sim,upx);
    auto &uyn = TF_getTmpField(sim,upx);
    auto &uzn = TF_getTmpField(sim,upx);

    tf2::oper_prod(M, ux, upx);
    tf2::oper_prod(M, uy, upy);
    tf2::oper_prod(M, uz, upz);

    tf2::oper_copy(upx,uxn);
    tf2::oper_copy(upy,uyn);
    tf2::oper_copy(upz,uzn);
    
    auto &diffx0 = tf2::getOrCreateField(sim,"diffx_0",upx);
    auto &diffy0 = tf2::getOrCreateField(sim,"diffy_0",upx);
    auto &diffz0 = tf2::getOrCreateField(sim,"diffz_0",upx);

    auto &convx0 = tf2::getOrCreateField(sim,"convx_0",upx);
    auto &convy0 = tf2::getOrCreateField(sim,"convy_0",upx);
    auto &convz0 = tf2::getOrCreateField(sim,"convz_0",upx);

    tf2::oper_prod(diffx0,diffx0,0.0);
    tf2::oper_prod(diffy0,diffy0,0.0);
    tf2::oper_prod(diffz0,diffz0,0.0);
    tf2::oper_prod(convx0,convx0,0.0);
    tf2::oper_prod(convy0,convy0,0.0);
    tf2::oper_prod(convz0,convz0,0.0);

    diffx.at(0) = &diffx0;
    diffy.at(0) = &diffy0;
    diffz.at(0) = &diffz0;

    convx.at(0) = &convx0;
    convy.at(0) = &convy0;
    convz.at(0) = &convz0; 

    //SUBSTAGE 1
    tf2::info("STAGE 1\n"); 
    // Predictor velocity should be computed here but for 1st stage it is not
    // required

    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
 
    diffx0  = diffusive(ux,0,sim);
    diffy0  = diffusive(uy,1,sim);
    diffz0  = diffusive(uz,2,sim);
    convx0  = convective(ufx,0,sim);
    convy0  = convective(ufy,1,sim);
    convz0  = convective(ufz,2,sim);

    for(uint64_t i=1; i<coefs.b.size(); ++i)
    {
        tf2::info("STAGE %d\n",i+1);
        auto &convxi = tf2::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
        auto &convyi = tf2::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
        auto &convzi = tf2::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
        auto &diffxi = tf2::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
        auto &diffyi = tf2::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
        auto &diffzi = tf2::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
    
        upx = predictorVector(uxn,convx,diffx,msx,coefs,dt,i,sim);
        upy = predictorVector(uyn,convy,diffy,msy,coefs,dt,i,sim);
        upz = predictorVector(uzn,convz,diffz,msz,coefs,dt,i,sim);
   
        tf2::oper_prod(ICF, upx, ufx);
        tf2::oper_prod(ICF, upy, ufy);
        tf2::oper_prod(ICF, upz, ufz);
    
        pSource = pRHS(ufx,ufy,ufz,sim); 
        tf2::oper_solve(psolver, pSource, p);
 
        projection(ufx,p,GX);
        projection(ufy,p,GY);
        projection(ufz,p,GZ);

        projection(upx,p,GX,IFC,sim);
        projection(upy,p,GY,IFC,sim);
        projection(upz,p,GZ,IFC,sim);
    
        // Map the velocity to the nodes.
        tf2::oper_prod(ID_CN, upx, ux);
        tf2::oper_prod(ID_CN, upy, uy);
        tf2::oper_prod(ID_CN, upz, uz);
 
        tf2::oper_prod(diffxi,diffxi,0.0);
        tf2::oper_prod(diffyi,diffyi,0.0);
        tf2::oper_prod(diffzi,diffzi,0.0);
        tf2::oper_prod(convxi,convxi,0.0);
        tf2::oper_prod(convyi,convyi,0.0);
        tf2::oper_prod(convzi,convzi,0.0);

        diffx.at(i) = &diffxi;
        diffy.at(i) = &diffyi;
        diffz.at(i) = &diffzi;
        convx.at(i) = &convxi;
        convy.at(i) = &convyi;
        convz.at(i) = &convzi;       

        diffxi  = diffusive(ux,0,sim);
        diffyi  = diffusive(uy,1,sim);
        diffzi  = diffusive(uz,2,sim);
        convxi  = convective(ufx,0,sim);
        convyi  = convective(ufy,1,sim);
        convzi  = convective(ufz,2,sim);
        
        tf2::info("Stage completed.\n");
    }

    upx = predictorVector(uxn,convx,diffx,msx,coefs,dt,0,sim);
    upy = predictorVector(uyn,convy,diffy,msy,coefs,dt,0,sim);
    upz = predictorVector(uzn,convz,diffz,msz,coefs,dt,0,sim);
   
    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pSource = pRHS(ufx,ufy,ufz,sim); 
    tf2::oper_solve(psolver, pSource, p);
 
    projection(ufx,p,GX);
    projection(ufy,p,GY);
    projection(ufz,p,GZ);

    projection(upx,p,GX,IFC,sim);
    projection(upy,p,GY,IFC,sim);
    projection(upz,p,GZ,IFC,sim);
    
    // Map the velocity to the nodes.
    tf2::oper_prod(ID_CN, upx, ux);
    tf2::oper_prod(ID_CN, upy, uy);
    tf2::oper_prod(ID_CN, upz, uz);
    tf2::info("Iteration completed.\n");
}

