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

TF_Func bool RK2Iteration_vector(tf2::Simulation &sim)
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

    std::vector<tf2::Field*> convx(2);
    std::vector<tf2::Field*> diffx(2);
    std::vector<tf2::Field*> convy(2);
    std::vector<tf2::Field*> diffy(2);
    std::vector<tf2::Field*> convz(2);
    std::vector<tf2::Field*> diffz(2);
    
    std::vector<double> b = {0.5,0.5};
    int s = b.size();
   
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

    diffx0  = diffusive(ux,0,sim);
    diffy0  = diffusive(uy,1,sim);
    diffz0  = diffusive(uz,2,sim);
    convx0  = convective(ufx,0,sim);
    convy0  = convective(ufy,1,sim);
    convz0  = convective(ufz,2,sim);

    //predictor 2 stage

    for(int i = 1; i<s; ++i){
      
      tf2::oper_copy(uxn,upx);
      tf2::oper_copy(uyn,upy);
      tf2::oper_copy(uzn,upz);

      for(int j=0; j<i; ++j){
        upx = predictor(upx,*(convx.at(j)),*(diffx.at(j)),msx,1.0,dt,sim);
        upy = predictor(upy,*(convy.at(j)),*(diffy.at(j)),msy,1.0,dt,sim);
        upz = predictor(upz,*(convz.at(j)),*(diffz.at(j)),msz,1.0,dt,sim);
      }
      
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
 
      auto &convxi = tf2::getOrCreateField(sim,"convx_"+std::to_string(i),upx);
      auto &convyi = tf2::getOrCreateField(sim,"convy_"+std::to_string(i),upx);
      auto &convzi = tf2::getOrCreateField(sim,"convz_"+std::to_string(i),upx);
      auto &diffxi = tf2::getOrCreateField(sim,"diffx_"+std::to_string(i),upx);
      auto &diffyi = tf2::getOrCreateField(sim,"diffy_"+std::to_string(i),upx);
      auto &diffzi = tf2::getOrCreateField(sim,"diffz_"+std::to_string(i),upx);
 
      diffx.at(1) = &diffxi;
      diffy.at(1) = &diffyi;
      diffz.at(1) = &diffzi;
      convx.at(1) = &convxi;
      convy.at(1) = &convyi;
      convz.at(1) = &convzi;       

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
    }
    //printMaxVals("end of first stage. DIFFUSIVE FIELD",diffx,diffy,diffz);
    //printMaxVals("end of first stage. CONVECTIVE FIELD",convx,convy,convz);
    
    //Predictor velocity
    tf2::oper_copy(uxn,upx);
    tf2::oper_copy(uyn,upy);
    tf2::oper_copy(uzn,upz);
    
    for(int i = 0; i < s; ++i){ 
      upx = predictor(upx,*(convx.at(i)),*(diffx.at(i)),msx,b.at(i),dt,sim);
      upy = predictor(upy,*(convy.at(i)),*(diffy.at(i)),msy,b.at(i),dt,sim);
      upz = predictor(upz,*(convz.at(i)),*(diffz.at(i)),msz,b.at(i),dt,sim);
    }

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
