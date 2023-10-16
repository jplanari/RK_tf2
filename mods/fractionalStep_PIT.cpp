#include "tf2/ll_Opers.h"
#include "tf2/Opers.h"
#include "tf2/Simulation.h"

TF_Func void SetUp_Momentum_PIT(tf2::Simulation &sim)
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
