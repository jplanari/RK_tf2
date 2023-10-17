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

TF_Func void SetUp_Momentum_RK(tf2::Simulation &sim)
{
    // This function exists basically to initialize all the fields required by
    // the official FSM, with the correct number of components, i.e., one
    // component per simulation.
    TF_uAssert((sim.meshes.size() == 1) xor tf2::hasSMesh(sim), "not yet implemented for multiple meshes");
    TF_uAssert(tf2::hasField(sim, "ux_N"), "required field ux_N not defined");
    TF_uAssert(tf2::hasField(sim, "uy_N"), "required field uy_N not defined");
    TF_uAssert(tf2::hasField(sim, "uz_N"), "required field uz_N not defined");
    TF_uAssert(tf2::hasField(sim, "Omega_C"), "required field Omega_C not defined");

    const auto &meshName = tf2::getMeshName(sim);
    const auto &msuf     = tf2::meshSuffix(meshName);
    const auto &cells    = "Cells"+msuf;
    const auto &faces    = "Faces"+msuf;
    const auto &nodes    = "Nodes"+msuf;

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

    const auto dim = tf2::getField(sim, "ux_N").dim;

    tf2::getOrCreateField(sim, dim, "upredx_C",  cells);
    tf2::getOrCreateField(sim, dim, "upredy_C",  cells);
    tf2::getOrCreateField(sim, dim, "upredz_C",  cells);
    tf2::getOrCreateField(sim, dim, "uxdiff_C",  cells);
    tf2::getOrCreateField(sim, dim, "uydiff_C",  cells);
    tf2::getOrCreateField(sim, dim, "uzdiff_C",  cells);
    tf2::getOrCreateField(sim, dim, "uxdiff0_C", cells);
    tf2::getOrCreateField(sim, dim, "uydiff0_C", cells);
    tf2::getOrCreateField(sim, dim, "uzdiff0_C", cells);
    tf2::getOrCreateField(sim, dim, "uxconv_C",  cells);
    tf2::getOrCreateField(sim, dim, "uyconv_C",  cells);
    tf2::getOrCreateField(sim, dim, "uzconv_C",  cells);
    tf2::getOrCreateField(sim, dim, "uxconv0_C", cells);
    tf2::getOrCreateField(sim, dim, "uyconv0_C", cells);
    tf2::getOrCreateField(sim, dim, "uzconv0_C", cells);
    tf2::getOrCreateField(sim, dim, "momSrcx_C", cells);
    tf2::getOrCreateField(sim, dim, "momSrcy_C", cells);
    tf2::getOrCreateField(sim, dim, "momSrcz_C", cells);
    tf2::getOrCreateField(sim, dim, "P_C",       cells);

    tf2::getOrCreateField(sim, dim, "ux_F", faces);
    tf2::getOrCreateField(sim, dim, "uy_F", faces);
    tf2::getOrCreateField(sim, dim, "uz_F", faces);

    tf2::getOrCreateField(sim, dim, "ux0_N", nodes);
    tf2::getOrCreateField(sim, dim, "uy0_N", nodes);
    tf2::getOrCreateField(sim, dim, "uz0_N", nodes);
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
 
    tf2::oper_prod(diffx0,diffx0,0.0);
    tf2::oper_prod(diffy0,diffy0,0.0);
    tf2::oper_prod(diffz0,diffz0,0.0);
    tf2::oper_prod(convx0,convx0,0.0);
    tf2::oper_prod(convy0,convy0,0.0);
    tf2::oper_prod(convz0,convz0,0.0);

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
    
        upx = predictorVector(uxn,convx,diffx,msx,coefs,dt,i+1,sim);
        upy = predictorVector(uyn,convy,diffy,msy,coefs,dt,i+1,sim);
        upz = predictorVector(uzn,convz,diffz,msz,coefs,dt,i+1,sim);
   
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
 
        diffx.at(i) = &diffxi;
        diffy.at(i) = &diffyi;
        diffz.at(i) = &diffzi;
        convx.at(i) = &convxi;
        convy.at(i) = &convyi;
        convz.at(i) = &convzi;       

        tf2::oper_prod(diffxi,diffxi,0.0);
        tf2::oper_prod(diffyi,diffyi,0.0);
        tf2::oper_prod(diffzi,diffzi,0.0);
        tf2::oper_prod(convxi,convxi,0.0);
        tf2::oper_prod(convyi,convyi,0.0);
        tf2::oper_prod(convzi,convzi,0.0);

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

