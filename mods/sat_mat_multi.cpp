#include "tf2/Opers.h"
#include "tf2/Simulation.h"

TF_Func void computeEV_real_mat(tf2::Simulation &sim)
{
    // NOTE: since this depends on the geometry only, it can be calculated once
    // at the start of the simulation.

    auto &ux_N = tf2::getField(sim, "ux_N");
    uint32_t dim = ux_N.dim;

    auto &vol = tf2::meshCellVolume(sim, dim);
    auto &af  = tf2::meshFaceArea(sim, dim);
    auto &dn  = tf2::meshFaceCellDist(sim, dim);
    auto &fim = tf2::meshFaceInnerMask(sim, dim);

    auto &SF  = tf2::getMatrix(sim, "SumFaces");

    auto &tmp = TF_getTmpField(sim, af);
    auto &gg  = TF_getTmpField(sim, SF, tmp);

    // tmp_F = (A/|d*n|)_F
    tf2::oper_copy(af, tmp);
    tf2::oper_div(tmp, dn);

    // If we do not want the boundary faces to contribute to the sum for each
    // cell, we just set the values of 'tmp' to zero at the boundary faces. We
    // achieve this by multiplying 'tmp' by 'fim' (face inner mask), where
    // 'fim' is 1 for inner faces (faces with exactly two neighbour cells) and
    // 0 for boundary faces.
    tf2::oper_prod(fim, tmp);

    // gg_C = (sum(faces) tmp_f)/vol_C
    tf2::oper_prod(SF, tmp, gg);
    tf2::oper_div(gg, vol);

    double kinVisc = 1.0/sim.IOParamD["Re"];
    sim.IOParamD["EVreal_mat"] = kinVisc*tf2::max(tf2::oper_max(gg));
}

TF_Func bool computeEV_imag_mat(tf2::Simulation &sim)
{
    // NOTE: since this depends on the velocity field, it has to be updated
    // every iteration.

    auto &ux_N = tf2::getField(sim, "ux_N");
    auto &uy_N = tf2::getField(sim, "uy_N");
    auto &uz_N = tf2::getField(sim, "uz_N");

    auto &INF = tf2::getMatrix(sim, "Interp_NF");

    auto &ux_f = TF_getTmpField(sim, INF, ux_N);
    auto &uy_f = TF_getTmpField(sim, INF, ux_N);
    auto &uz_f = TF_getTmpField(sim, INF, ux_N);

    tf2::oper_prod(INF, ux_N, ux_f);
    tf2::oper_prod(INF, uy_N, uy_f);
    tf2::oper_prod(INF, uz_N, uz_f);

    // In a real simulation using the fractional step method, we will have the
    // velocity at the faces already, so we just get it. There is no need to
    // interpolate the velocity from the nodes.
    // auto &ux_f = tf2::getField(sim, "ux_F");
    // auto &uy_f = tf2::getField(sim, "uy_F");
    // auto &uz_f = tf2::getField(sim, "uz_F");

    uint32_t dim = ux_N.dim;

    auto &nx  = tf2::meshFaceNx(sim, dim);
    auto &ny  = tf2::meshFaceNy(sim, dim);
    auto &nz  = tf2::meshFaceNz(sim, dim);
    auto &af  = tf2::meshFaceArea(sim, dim);
    auto &fim = tf2::meshFaceInnerMask(sim, dim);
    auto &vol = tf2::meshCellVolume(sim, dim);

    auto &SF  = tf2::getMatrix(sim, "SumFaces");
    auto &tmp = TF_getTmpField(sim, ux_f);
    auto &gg  = TF_getTmpField(sim, SF, tmp);

    // tmp_F = (A*|u·n|)_F
    tf2::oper_fmadd(ux_f, nx, tmp, 1.0, 0.0);
    tf2::oper_fmadd(uy_f, ny, tmp, 1.0, 1.0);
    tf2::oper_fmadd(uz_f, nz, tmp, 1.0, 1.0);
    tf2::oper_max(tmp, tmp, 1.0, -1.0);
    tf2::oper_prod(af, tmp);

    // If we do not want the boundary faces to contribute to the sum for each
    // cell, we just set the values of 'tmp' to zero at the boundary faces. We
    // achieve this by multiplying 'tmp' by 'fim' (face inner mask), where
    // 'fim' is 1 for inner faces (faces with exactly two neighbour cells) and
    // 0 for boundary faces.
    tf2::oper_prod(fim, tmp);

    // gg_C = (sum(faces) tmp_f)/vol_C
    tf2::oper_prod(SF, tmp, gg);
    tf2::oper_div(gg, vol);

    sim.IOParamD["EVimag_mat"] = 0.5*tf2::max(tf2::oper_max(gg));
    return tf2::Iter_Continue;
}

TF_Func void compare_ev(tf2::Simulation &sim)
{
    tf2::info("evR:    %.10e\n", sim.IOParamD["_EVreal"]);
    tf2::info("evRmat: %.10e\n", sim.IOParamD["EVreal_mat"]);
    tf2::info("evI:    %.10e\n", sim.IOParamD["_EVimag"]);
    tf2::info("evImat: %.10e\n", sim.IOParamD["EVimag_mat"]);
}

