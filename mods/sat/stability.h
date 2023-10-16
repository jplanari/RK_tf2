double tauMethod(int s, double f)
{
    return 1 + (s-1)*(1-f);
}

double stabCoefs(int i, butcherTableau coefs)
{
    double sum=0.0;
    double isum1=0.0;
    double isum2=0.0;
    double isum3=0.0;
    double isum4=0.0;
    switch(i)
    {
	case 0:
	    sum = 1.0;
	    break;
	case 1:
	    for(int j=0;j<(int)coefs.b.size();++j) sum += coefs.b.at(j);
	    break;
	case 2:
	    for(int j=1;j<(int)coefs.b.size();++j){
	        isum1=0.0;
		for(int k=1;k<=j;++k) isum1 += coefs.A.at(id(j+1,k));
		sum += isum1*coefs.b.at(j);
	    }
	    break;
	case 3:
	    for(int j=2;j<(int)coefs.b.size();++j){
	        isum1=0.0;
		for(int k=2;k<=j;++k){
		     isum2=0.0;
		     for(int l=1;l<k;++l) isum2 += coefs.A.at(id(k,l));
		     isum1 += isum2*coefs.A.at(id(j+1,k));
		}
		sum += isum1*coefs.b.at(j);
	    }
	    break;
	case 4:
	    for(int j=3;j<(int)coefs.b.size();++j){
		isum1=0.0;
		for(int k=3;k<=j;++k){
		    isum2=0.0;
		    for(int l=2;l<k;++l){
			isum3=0.0;
			for(int m=1;m<l;++m) isum3 += coefs.A.at(id(l,m));
			isum2 += isum3*coefs.A.at(id(k,l));
		    }
		    isum1 += isum2*coefs.A.at(id(j+1,k));
		}
		sum += isum1*coefs.b.at(j);
	    }
	    break;
	case 5:
	    for(int j=4;j<(int)coefs.b.size();++j){
		isum1=0.0;
		for(int k=4;k<=j;++k){
		    isum2=0.0;
		    for(int l=3;l<k;++l){
			isum3 = 0.0;
			for(int m=2; m<l; ++m){
			    isum4 = 0.0;
			    for(int n=1; n<m; ++n) isum4 += coefs.A.at(id(m,n));
			    isum3 += isum4*coefs.A.at(id(l,m));
			}
			isum2 += isum3*coefs.A.at(id(k,l));
		    }
		    isum1 += isum2*coefs.A.at(id(j+1,k));
		}
		sum += isum1*coefs.b.at(j);
	    }
	    break;
    }
    return sum;
}

double stabilityRK_coefs(double h, double phi, butcherTableau coefs)
{
    double real=0.0, imag=0.0, coef;
    for(int i=0; i<=(int)coefs.b.size(); ++i)
    {
	coef = stabCoefs(i,coefs);
	real += coef*pow(h,i)*cos(i*phi);
        imag += coef*pow(h,i)*sin(i*phi);
    }
    return sqrt(real*real+imag*imag);
}

double stabilityRegion(double phi, butcherTableau coefs,  double param=0.0)
{
    double hmax = 6;
    double h=0.01;
    double dh=0.01;
    double ra=0.0,prevra=2.0;
    while (h<hmax)
    {
        ra = stabilityRK_coefs(h,phi,coefs);
        if (ra>1.0 && prevra<1.0) break;
        else prevra=ra;
        h += dh;
    }
    if (h<(hmax-dh)) return h;
    else return 0.0;
}


