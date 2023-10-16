#include "tf2/ll_Opers.h"
#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include "fsm_rk/butcherTableaux.h"
#include "fsm_rk/setup_fsm.h"
#include "fsm_rk/cd.h"
#include "fsm_rk/fsm_rk.h"
#include "fsm_rk/hardcoded.h"

#include <vector>
#include <cstdarg>

// +---------------------------------------------------------------------------+
// |                              General remarks                              |
// +---------------------------------------------------------------------------+

// This module introduces some functions to help solving the Navier-Stokes
// equations. For the moment only incompressible fluids have been considered.


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
    
    butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

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

