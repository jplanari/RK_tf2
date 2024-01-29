#include "tf2/ll_Opers.h"
#include "tf2/Opers.h"
#include "tf2/Simulation.h"

#include "fsm_rk/butcherTableaux.h"
#include "fsm_rk/setup_fsm.h"
#include "fsm_rk/cd.h"
#include "fsm_rk/fsm_rk.h"

#include "energy/energy.h"

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
    if(sim.IOParamD["_ElapsedTime"]>sim.IOParamD["_MaxTime"]) return tf2::Iter_Stop;
    
    auto &M     = tf2::getMatrix(sim, "Id_NC");
    auto &GX    = tf2::getMatrix(sim, "GradX_CF");
    auto &GY    = tf2::getMatrix(sim, "GradY_CF");
    auto &GZ    = tf2::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tf2::getMatrix(sim, "Interp_CF");
    auto &IFC   = tf2::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
    auto &INF   = tf2::getMatrix(sim, "Interp_NF");

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

    std::vector<double> b = coefs.b;
    std::vector<double> A = coefs.A;
    int s = b.size();
 
    std::vector<tf2::Field*> convx(s);
    std::vector<tf2::Field*> diffx(s);
    std::vector<tf2::Field*> convy(s);
    std::vector<tf2::Field*> diffy(s);
    std::vector<tf2::Field*> convz(s);
    std::vector<tf2::Field*> diffz(s);  

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

    diffusive(ux,diffx0,sim);
    diffusive(uy,diffy0,sim);
    diffusive(uz,diffz0,sim);
    convective(ufx,convx0,sim);
    convective(ufy,convy0,sim);
    convective(ufz,convz0,sim);

    //predictor 2 stage

    for(int i = 1; i<s; ++i){
      
      tf2::oper_copy(uxn,upx);
      tf2::oper_copy(uyn,upy);
      tf2::oper_copy(uzn,upz);

      for(int j=0; j<i; ++j){
        predictor(upx,upx,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt,sim);
        predictor(upy,upy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt,sim);
        predictor(upz,upz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt,sim);
      }
      
      tf2::oper_prod(ICF, upx, ufx);
      tf2::oper_prod(ICF, upy, ufy);
      tf2::oper_prod(ICF, upz, ufz);
    
      pRHS(pSource,ufx,ufy,ufz,sim); 
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
 
      auto &convxi = tf2::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
      auto &convyi = tf2::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
      auto &convzi = tf2::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
      auto &diffxi = tf2::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
      auto &diffyi = tf2::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
      auto &diffzi = tf2::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
 
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

      diffusive(ux,diffxi,sim);
      diffusive(uy,diffyi,sim);
      diffusive(uz,diffzi,sim);
      convective(ufx,convxi,sim);
      convective(ufy,convyi,sim);
      convective(ufz,convzi,sim);
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tf2::oper_copy(uxn,upx);
    tf2::oper_copy(uyn,upy);
    tf2::oper_copy(uzn,upz);
    
    for(int i = 0; i < s; ++i){ 
      predictor(upx,upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt,sim);
      predictor(upy,upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt,sim);
      predictor(upz,upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt,sim);
    }

    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pRHS(pSource,ufx,ufy,ufz,sim); 
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

    return tf2::Iter_Continue;
}


TF_Func bool RKiteration_energy(tf2::Simulation &sim)
{
    if(sim.IOParamD["_ElapsedTime"]>sim.IOParamD["_MaxTime"]) return tf2::Iter_Stop;
    
    auto &M     = tf2::getMatrix(sim, "Id_NC");
    auto &GX    = tf2::getMatrix(sim, "GradX_CF");
    auto &GY    = tf2::getMatrix(sim, "GradY_CF");
    auto &GZ    = tf2::getMatrix(sim, "GradZ_CF");
    auto &ICF   = tf2::getMatrix(sim, "Interp_CF");
    auto &IFC   = tf2::getMatrix(sim, "Interp_FC");
    auto &ID_CN = tf2::getMatrix(sim, "Id_CN");
    auto &INF   = tf2::getMatrix(sim, "Interp_NF");
    
    auto &ux     = tf2::getField(sim, "ux_N");
    auto &uy     = tf2::getField(sim, "uy_N");
    auto &uz     = tf2::getField(sim, "uz_N");
    
    auto &T      = tf2::getField(sim, "T_N");
    
    auto &ufx    = tf2::getField(sim, "ux_F");
    auto &ufy    = tf2::getField(sim, "uy_F");
    auto &ufz    = tf2::getField(sim, "uz_F");

    auto &upx    = tf2::getField(sim, "upredx_C");
    auto &upy    = tf2::getField(sim, "upredy_C");
    auto &upz    = tf2::getField(sim, "upredz_C");

    auto &msx    = tf2::getField(sim, "momSrcx_C");
    auto &msy    = tf2::getField(sim, "momSrcy_C");
    auto &msz    = tf2::getField(sim, "momSrcz_C");
    auto &msT    = TF_getTmpField(sim, msx);

    tf2::oper_setConst(msT,0.0);
    
    auto &p      = tf2::getField(sim, "P_C");
    static auto psolver = tf2::getSolver(sim, "Pressure_Solver");
    butcherTableau coefs = intScheme(sim.IOParamS["RKmethod"]);

    double lambda = sim.IOParamD["lambda"]; //change when set properly

    std::vector<double> b = coefs.b;
    std::vector<double> A = coefs.A;
    int s = b.size();
 
    std::vector<tf2::Field*> convx(s);
    std::vector<tf2::Field*> diffx(s);
    std::vector<tf2::Field*> convy(s);
    std::vector<tf2::Field*> diffy(s);
    std::vector<tf2::Field*> convz(s);
    std::vector<tf2::Field*> diffz(s);  

    std::vector<tf2::Field*> Trk(s); //stores the temperature at the cells for the Runge-Kutta substages
    
    // Temporary variables.
    auto &pSource = TF_getTmpField(sim, p);
    
    auto &uxn = TF_getTmpField(sim,upx);
    auto &uyn = TF_getTmpField(sim,upx);
    auto &uzn = TF_getTmpField(sim,upx);

    auto &Tn = TF_getTmpField(sim,T);
    auto &Tc = TF_getTmpField(sim,upx);
    auto &Tf = TF_getTmpField(sim,ufx);
    
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
   
    auto &diffx0 = tf2::getOrCreateField(sim,"diffx_0",upx);
    auto &diffy0 = tf2::getOrCreateField(sim,"diffy_0",upx);
    auto &diffz0 = tf2::getOrCreateField(sim,"diffz_0",upx);

    auto &convx0 = tf2::getOrCreateField(sim,"convx_0",upx);
    auto &convy0 = tf2::getOrCreateField(sim,"convy_0",upx);
    auto &convz0 = tf2::getOrCreateField(sim,"convz_0",upx);

    auto &diffT = tf2::getOrCreateField(sim,"diffT",upx);
    auto &convT = tf2::getOrCreateField(sim,"convT",upx);
    auto &T0 = tf2::getOrCreateField(sim,"T_0",T);

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

    diffusive(ux,diffx0,sim);
    diffusive(uy,diffy0,sim);
    diffusive(uz,diffz0,sim);
    convective(ufx,convx0,sim);
    convective(ufy,convy0,sim);
    convective(ufz,convz0,sim);

    tf2::oper_prod(T0,T0,0.0);
    Trk.at(0) = &T0;
    tf2::oper_copy(T,T0);
    
    //predictor 2 stage

    for(int i = 1; i<s; ++i){
      
      tf2::oper_copy(uxn,upx);
      tf2::oper_copy(uyn,upy);
      tf2::oper_copy(uzn,upz);

      for(int j=0; j<i; ++j){
        generateSourceTerm(0.0,*(Trk.at(j)),msx,sim); 
        predictor(upx,upx,*(convx.at(j)),*(diffx.at(j)),msx,A.at(id(i+1,j+1)),dt,sim);
        generateSourceTerm(1.0,*(Trk.at(j)),msy,sim); 
        predictor(upy,upy,*(convy.at(j)),*(diffy.at(j)),msy,A.at(id(i+1,j+1)),dt,sim);
        generateSourceTerm(0.0,*(Trk.at(j)),msz,sim); 
        predictor(upz,upz,*(convz.at(j)),*(diffz.at(j)),msz,A.at(id(i+1,j+1)),dt,sim);
      }
      
      tf2::oper_prod(ICF, upx, ufx);
      tf2::oper_prod(ICF, upy, ufy);
      tf2::oper_prod(ICF, upz, ufz);
    
      pRHS(pSource,ufx,ufy,ufz,sim); 
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

      //ENERGY EQUATION
      
      tf2::oper_prod(M,T0,Tc);
      
      for(int j=0; j<i; ++j){
        tf2::oper_prod(INF,*(Trk.at(j)),Tf);
        convective(Tf,convT,sim);
        diffusive(*(Trk.at(j)),diffT,lambda,sim);
        predictor(Tc,Tc,convT,diffT,msT,A.at(id(i+1,j+1)),dt,sim);
      }
      
      auto &Ti = tf2::getOrCreateField(sim,"T_"+std::to_string(i),T0);

      Trk.at(i) = &Ti;
      tf2::oper_prod(Ti,Ti,0.0);
      tf2::oper_prod(ID_CN, Tc, Ti);

      auto &convxi = tf2::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
      auto &convyi = tf2::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
      auto &convzi = tf2::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
      auto &diffxi = tf2::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
      auto &diffyi = tf2::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
      auto &diffzi = tf2::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
 
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

      diffusive(ux,diffxi,sim);
      diffusive(uy,diffyi,sim);
      diffusive(uz,diffzi,sim);
      convective(ufx,convxi,sim);
      convective(ufy,convyi,sim);
      convective(ufz,convzi,sim);
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tf2::oper_copy(uxn,upx);
    tf2::oper_copy(uyn,upy);
    tf2::oper_copy(uzn,upz);
    
    for(int i = 0; i < s; ++i){ 
      generateSourceTerm(0.0,*(Trk.at(i)),msx,sim); 
      predictor(upx,upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt,sim);
      generateSourceTerm(1.0,*(Trk.at(i)),msy,sim); 
      predictor(upy,upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt,sim);
      generateSourceTerm(0.0,*(Trk.at(i)),msz,sim); 
      predictor(upz,upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt,sim);
    }

    tf2::oper_prod(ICF, upx, ufx);
    tf2::oper_prod(ICF, upy, ufy);
    tf2::oper_prod(ICF, upz, ufz);
    
    pRHS(pSource,ufx,ufy,ufz,sim); 
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

    for(int i = 0; i < s; ++i){ 
      tf2::oper_prod(INF,*(Trk.at(i)),Tf);
      convective(Tf,convT,sim);
      diffusive(*(Trk.at(i)),diffT,lambda,sim);
      predictor(Tc,Tc,convT,diffT,msT,b.at(i),dt,sim);
    }

    tf2::oper_prod(ID_CN, Tc, T);

    return tf2::Iter_Continue;
    //printMaxVals("end of time-step",ux,uy,uz);
}


