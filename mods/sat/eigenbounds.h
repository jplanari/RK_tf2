#include "meshFcs.h"

double computeUFaceContributionGershgorin_real(tf2::UMesh &m, mt::FaceIdT f, mt::CellIdT c)
{
    double real=0.0;
    auto cntr = m.getCellCentroid(c);
    auto c_nb = m.getOtherCell(f,c);
    if (c_nb != mt::CellIdT::None)
    {
        auto cntr_nb = m.getCellCentroid(c_nb);
        tf2::V3 dvec = cntr-cntr_nb;
        real = fabs(m.getFaceArea(f)/(dvec*m.getFaceNormal(f,c)));
    }
    return real;
}

double computeSFaceContributionGershgorin_real(tf2::SMesh &m, tf2::IJK F, tf2::IJK C, tf2::IJK CNB)
{
    auto cntr = tf2::calcCellCentroid(m,C);
    auto cntr_nb = tf2::calcCellCentroid(m,CNB);
    tf2::V3 dvec = cntr-cntr_nb;
    auto fnv = tf2::calcFaceNormal(m,F);
    double area = fnv*tf2::calcCellFacesArea(m,C);

    return fabs(area/(dvec*fnv));
}

double computeFaceContributionGershgorin_imag(tf2::V3 fnv, tf2::V3 uf, double Af)
{
    return 0.5*fabs(fnv*uf)*Af;
}

void computeImagEV_Gershgorin(tf2::Simulation &sim)
{
    auto &ux_N = tf2::getField(sim,"ux_N");
    auto &uy_N = tf2::getField(sim,"uy_N");
    auto &uz_N = tf2::getField(sim,"uz_N");

    auto &INF = tf2::getMatrix(sim,"Interp_NF");
    
    auto &uxf = TF_getTmpField(sim,ux_N.dim,"Faces");
    auto &uyf = TF_getTmpField(sim,ux_N.dim,"Faces");
    auto &uzf = TF_getTmpField(sim,ux_N.dim,"Faces");

    tf2::oper_prod(INF,ux_N,uxf);
    tf2::oper_prod(INF,uy_N,uyf);
    tf2::oper_prod(INF,uz_N,uzf);

    double u,v,w;
    tf2::V3 uf;
    double gershContImag;
    double evi=0.0;
    if(!tf2::hasSMesh(sim)){
        tf2::UMesh &m = tf2::getUMesh(sim);
    
    	  for (auto c : m.getCellsRange()){
          gershContImag = 0.0;

        	for (auto f : m.getCellFaces(c)){
            auto fnv = m.getFaceNormal(f, c);

	          u = uxf.array[TO_LOC(f)];
        	  v = uyf.array[TO_LOC(f)];
	          w = uzf.array[TO_LOC(f)];

        	  uf = {u,v,w};

            gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,m.getFaceArea(f))/m.getCellVolume(c);
        	}
        	
		      if (gershContImag>=evi) evi = gershContImag;
    	  }

    } else {
	    tf2::SMesh &m = tf2::getSMesh(sim);
      tf2::IJK F,C;
	    auto offset = m.topo.firstOwnedFaces[tf2::mpiRank()];
      uint32_t f;
	    tf2::V3 cfa,fnv;
	    for(uint32_t c = 0; c < tf2::getNumCells(m); ++c){
		
		    gershContImag = 0.0;
		    C = tf2::getCellIJKFromId(m,c);
		    F = C;
		
	      cfa = tf2::calcCellFacesArea(m,C);
		
		    for(uint32_t fc = 0; fc < 6; fc++){
			  
          if(hasNBCell(m,C,fc)){
            F.f = fc;
				    f = tf2::getFaceIndex(m,C.i,C.j,C.k,fc)-offset;
				    fnv = tf2::calcFaceNormal(m,F);
			
				    u = uxf.array[TO_LOC(f)];
				    v = uyf.array[TO_LOC(f)];
				    w = uzf.array[TO_LOC(f)];

				    uf = {u,v,w};
				    gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,fnv*cfa)/tf2::calcCellVolume(m,C);			
			    }
		    }
		    evi = std::max(evi,gershContImag);
	    }
    }
    sim.IOParamD["_EVimag"] = tf2::allReduce(evi, MPI_MAX);
}

TF_Func void computeEV_Gershgorin(tf2::Simulation &sim)
{
    auto &ux_N = tf2::getField(sim,"ux_N");
    auto &uy_N = tf2::getField(sim,"uy_N");
    auto &uz_N = tf2::getField(sim,"uz_N");
    auto &INF = tf2::getMatrix(sim,"Interp_NF");
    
    auto &ux_f = TF_getTmpField(sim,ux_N.dim,"Faces");
    auto &uy_f = TF_getTmpField(sim,ux_N.dim,"Faces");
    auto &uz_f = TF_getTmpField(sim,ux_N.dim,"Faces");

    tf2::oper_prod(INF,ux_N,ux_f);
    tf2::oper_prod(INF,uy_N,uy_f);
    tf2::oper_prod(INF,uz_N,uz_f);

    double kinVisc = sim.IOParamD["kinVisc"];

    double gershContImag;
    double gershContReal;

    double evr=0.0,evi=0.0;

    double u,v,w;
    tf2::V3 uf;

    if (!hasSMesh(sim)){

      tf2::UMesh &m = tf2::getUMesh(sim);
      for (auto c : m.getCellsRange()){
        
        gershContImag = 0.0;
        gershContReal = 0.0;

        for (auto f : m.getCellFaces(c)){

          auto fnv = m.getFaceNormal(f, c);

          u = ux_f.array[TO_LOC(f)];
          v = uy_f.array[TO_LOC(f)];
          w = uz_f.array[TO_LOC(f)];

          uf = {u,v,w};

          gershContReal += 2.0*kinVisc*computeUFaceContributionGershgorin_real(m,f,c)/m.getCellVolume(c);
          gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,m.getFaceArea(f))/m.getCellVolume(c);
        }
	      evr = std::max(evr,gershContReal);
	      evi = std::max(evi,gershContImag);
      }
    } else {
	    tf2::SMesh &m = tf2::getSMesh(sim);
      tf2::IJK F,C,CNB;
	    auto offset = m.topo.firstOwnedFaces[tf2::mpiRank()];
      uint32_t f;
	    tf2::V3 fnv, cfa;
	    for(uint32_t c = 0; c < tf2::getNumCells(m); ++c) {
		    gershContImag = 0.0;
		    gershContReal = 0.0;
		
		    C = tf2::getCellIJKFromId(m,c);
		    F = C;
		    CNB = tf2::getCellIJKFromId(m,c);
		    cfa = tf2::calcCellFacesArea(m,C);
		    for(uint32_t fc = 0; fc < 6; fc++){
			    if(hasNBCell(m,C,fc)){
				    f = tf2::getFaceIndex(m,C.i,C.j,C.k,fc)-offset;
				    F.f = fc;
				    CNB = C;
				    fnv = tf2::calcFaceNormal(m,F);
				
				    switch(fc){
				      case 0:
				    	  CNB.i = C.i-1;
					      break;
				      case 1:
				    	  CNB.i = C.i+1;
					      break;
				      case 2:
				        CNB.j = C.j-1;
					      break;
				      case 3:
				        CNB.j = C.j+1;
					      break;
				      case 4:
				        CNB.k = C.k-1;
					      break;
				      case 5:
					      CNB.k = C.k+1;
					      break;
				    }
			
				    u = ux_f.array[TO_LOC(f)];
				    v = uy_f.array[TO_LOC(f)];
				    w = uz_f.array[TO_LOC(f)];

				    uf = {u,v,w};
				    gershContReal += 2.0*kinVisc*computeSFaceContributionGershgorin_real(m,F,C,CNB)/tf2::calcCellVolume(m,C);
				    gershContImag += computeFaceContributionGershgorin_imag(fnv,uf,fnv*cfa)/tf2::calcCellVolume(m,C);			
			    }
		    }
        evr = std::max(evr,gershContReal);
		    evi = std::max(evi,gershContImag);
	    }
    }
    sim.IOParamD["_EVreal"] = tf2::allReduce(evr, MPI_MAX);
    sim.IOParamD["_EVimag"] = tf2::allReduce(evi, MPI_MAX);
}

TF_Func void computeRealEV_GershgorinMat(tf2::Simulation &sim)
{
    // NOTE: since this depends on the geometry only, it can be calculated once
    // at the start of the simulation.

    auto &ux_N = tf2::getField(sim, "ux_N");
    uint32_t dim = ux_N.dim;

    auto &vol = tf2::meshCellVolume(sim, dim);
    auto &af  = tf2::meshFaceArea(sim, dim);
    auto &dn  = faceCellDist(sim);
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

    double kinVisc = sim.IOParamD["kinVisc"];
    sim.IOParamD["_EVreal"] = kinVisc*tf2::max(tf2::oper_max(gg));
}

TF_Func bool computeImagEV_GershgorinMat(tf2::Simulation &sim)
{
    // NOTE: since this depends on the velocity field, it has to be updated
    // every iteration.

    auto &ux_f = tf2::getField(sim, "ux_F");
    auto &uy_f = tf2::getField(sim, "uy_F");
    auto &uz_f = tf2::getField(sim, "uz_F");

    // In a real simulation using the fractional step method, we will have the
    // velocity at the faces already, so we just get it. There is no need to
    // interpolate the velocity from the nodes.
    // auto &ux_f = tf2::getField(sim, "ux_F");
    // auto &uy_f = tf2::getField(sim, "uy_F");
    // auto &uz_f = tf2::getField(sim, "uz_F");

    uint32_t dim = ux_f.dim;

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

    sim.IOParamD["_EVimag"] = 0.5*tf2::max(tf2::oper_max(gg));
    return tf2::Iter_Continue;
}

tf2::Field &computeLambdaTilde(tf2::Simulation &sim)
{
    auto dim = tf2::getField(sim,"ux_N").dim;
    auto &result = tf2::getOrCreateField(sim, dim, "Lambda", "Faces");
    
    auto &af = tf2::meshFaceArea(sim,dim);
    auto &vol = calcFaceStaggeredVolumes(sim);
    

    tf2::oper_axpy(af,result,1.0,0.0);
    tf2::oper_prod(af,result,sim.IOParamD["kinVisc"]);
    tf2::oper_div(result,vol);
    
    return result;
}

tf2::Field &computeAbsMassFluxes(tf2::Simulation &sim)
{
    auto &ux_f = tf2::getField(sim, "ux_F");
    auto &uy_f = tf2::getField(sim, "uy_F");
    auto &uz_f = tf2::getField(sim, "uz_F");
    auto dim = ux_f.dim;

    auto &nx  = tf2::meshFaceNx(sim,dim);
    auto &ny  = tf2::meshFaceNy(sim,dim);
    auto &nz  = tf2::meshFaceNz(sim,dim);
    auto &af  = tf2::meshFaceArea(sim,dim);

    auto &result = tf2::getOrCreateField(sim, dim, "Fs", "Faces");
    
    tf2::oper_fmadd(ux_f, nx, result, 1.0, 0.0);
    tf2::oper_fmadd(uy_f, ny, result, 1.0, 1.0);
    tf2::oper_fmadd(uz_f, nz, result, 1.0, 1.0);
    tf2::oper_prod(af,result);
    tf2::oper_apply(result,[] (double x) {return fabs(x);});

    return result;    
}

double compute_AlgEigCD(tf2::Field &Fs, double factor, tf2::Simulation &sim)
{
    uint32_t resultSize = tf2::getNumEntries(Fs);
    double eigenbound = 0.0; 
    if(tf2::hasSMesh(sim))
    {
      tf2::IJK F,C;
      const auto &m = tf2::getSMesh(sim);
	    auto offset = m.topo.firstOwnedFaces[tf2::mpiRank()];
      double contr, c_contr;
      for (uint32_t it = 0; it < resultSize; it++)
      {
        F = tf2::getFaceIJKFromId(m,it);
        std::vector<tf2::IJK> cnb = getFaceNbCells(F,m,sim);
        tf2::IJK CNB;
        contr=0.0;
        for(int c=0; c<cnb.size(); c++)
        {
          CNB = cnb.at(c);
          c_contr=0.0;
          for(int f=0; f<6; ++f)
          {
            auto ff = tf2::getFaceIndex(m,CNB.i,CNB.j,CNB.k,f)-offset;
            auto FF = tf2::getFaceIJKFromId(m,ff);
            if(!tf2::isBoundaryFace(m,FF))
              c_contr += Fs.array[TO_U32(ff)];
          }
          contr += c_contr/tf2::calcCellVolume(m,CNB);
        }
        eigenbound = std::max(eigenbound,contr);
      }
    }
    return factor*eigenbound;    
}

TF_Func bool computeImagEV_AlgEigCD(tf2::Simulation &sim)
{
  auto &Fs = computeAbsMassFluxes(sim);
  sim.IOParamD["_EVimag"] = compute_AlgEigCD(Fs,0.25,sim);
  return tf2::Iter_Continue;
}

TF_Func void computeRealEV_AlgEigCD(tf2::Simulation &sim)
{
  auto &Fs = computeLambdaTilde(sim);
  sim.IOParamD["_EVreal"] = compute_AlgEigCD(Fs,1.0,sim);
}

TF_Func void computeEV_AlgEigCD(tf2::Simulation &sim)
{
  computeRealEV_AlgEigCD(sim);
  bool imag = computeImagEV_AlgEigCD(sim);
}

TF_Func void printEV(tf2::Simulation &sim)
{
  tf2::info("Eigenvalues:\n");
  tf2::info("\tReal = %e\n",sim.IOParamD["_EVreal"]);
  tf2::info("\tImag = %e\n",sim.IOParamD["_EVimag"]);
}
