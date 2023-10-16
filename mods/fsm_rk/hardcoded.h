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
