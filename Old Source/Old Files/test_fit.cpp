#include "simulator.h"

int main(int argc,char** argv){

  string outfile("output.fits");
  string modfile("model.fits");
  string obsfile("observation.fits");
  
  if(argc > 1)
    outfile = argv[1];

  int pnum,znum;
  FITS *pInfile,*pInfile2;
  pInfile = new FITS(modfile,Read);
  pInfile2 = new FITS(obsfile,Read);
  
  //reading primary header from fits file
  HDU& table = pInfile->pHDU();
  HDU& table2 = pInfile2->pHDU();
  //get number of parameters and initialize dynamic memory
  table.readKey("P_NUM",pnum);
  table.readKey("Z_NUM",znum);

  string pi[] = {"0","1","2","3","4","5"};

  double zmean[znum],zsigma[znum],zrange[2];
  
  printf("\n%s\n","Initial Redshift Parameters:");
  for(int i=0;i<znum;i++){
    table.readKey("Z"+pi[i]+"_MEAN",zmean[i]);
    table.readKey("Z"+pi[i]+"_SIGMA",zsigma[i]);
    printf("Mean: %5.3f, Sigma: %5.3f \n",zmean[i],zsigma[i]);
  }
  table.readKey("Z0_MIN",zrange[0]);
  table.readKey("Z0_MAX",zrange[1]);
  
  double bs[3],errs[3],flims[3];

  for(int i=0;i<3;i++){
    table2.readKey("WAVE_"+pi[i+1],bs[i]);
    table2.readKey("W"+pi[i+1]+"_FMIN",flims[i]);
    table2.readKey("W"+pi[i+1]+"_FERR",errs[i]);
    flims[i]*=1e-3;
    errs[i]*=1e-3;
  }

  //Make generic to any three bands
  printf("\n%s\n","Luminosity Function Parameters");
  lumfunct lums(flims[1]);

  double temp;
  table.readKey("PHI0",temp);
  printf("%8s: %6.3f\n","Phi0",pow(10,temp));
  lums.set_phi0(pow(10,temp));
  table.readKey("L0",temp);
  printf("%8s: %6.3f\n","L0",temp);
  lums.set_L0(temp);
  table.readKey("ALPHA",temp);
  printf("%8s: %6.3f\n","Alpha",temp);
  lums.set_alpha(temp);
  table.readKey("BETA",temp);
  printf("%8s: %6.3f\n","Beta",temp);
  lums.set_beta(temp);
  table.readKey("P",temp);
  printf("%8s: %6.3f\n","P",temp);
  lums.set_p(temp);
  table.readKey("Q",temp);
  printf("%8s: %6.3f\n","Q",temp);
  lums.set_q(temp);
  table.readKey("LAMBDA0",temp);
  printf("%8s: %5.3f\n","Lambda0",temp);
  lums.set_Lambda0(temp);
  table.readKey("OMEGAM",temp);
  printf("%8s: %5.3f\n","OmegaM",temp);
  lums.set_OmegaM(temp);
  table.readKey("H0",temp);
  printf("%8s: %7.3f\n","H0",temp);
  lums.set_H0(temp);
  
  lums.initialize();

  simulator tester(bs,errs,flims,&lums,modfile,obsfile);
  tester.set_zrange(zrange);

  printf("\n%s\n","Initial Model Parameters:");
  fixed_params ps;
  ps.pnum = pnum;
  ps.znum = znum;
  ps.p = new mparam[pnum];
  ps.sim = &tester;
  for (int i=0;i<pnum;i++){
    ps.p[i].p_id = i;
    table.readKey("P"+pi[i]+"_MEAN",ps.p[i].mean);
    table.readKey("P"+pi[i]+"_SIGMA",ps.p[i].sigma);
    table.readKey("P"+pi[i]+"_MIN",ps.p[i].min);
    table.readKey("P"+pi[i]+"_MAX",ps.p[i].max);
    printf("Mean: %5.3f, Sigma: %5.3f\n",ps.p[i].mean,ps.p[i].sigma);
  }

  delete pInfile;
  delete pInfile2;

  const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
  gsl_multimin_fminimizer *s = NULL;
  gsl_multimin_function minex_func;
  
  size_t iter = 0;
  int status;
  double size;

  int vecsize = 2;
  if(znum > 1)
    vecsize = 5;

  double chi_min = 1000;
  gsl_vector *best= gsl_vector_alloc(vecsize);
  gsl_vector *x= gsl_vector_alloc(vecsize);
  gsl_vector *ss= gsl_vector_alloc(vecsize);
  
  minex_func.n = vecsize;
  minex_func.f = simulate;
  minex_func.params = &ps;

  printf("\n%s\n","_______________________________________________________________");

  for (int zi=0;zi<znum;zi++){
    gsl_vector_set(best,0+zi*3,zmean[zi]);
    gsl_vector_set(best,1+zi*3,zsigma[zi]);
    if(zi > 0)
      gsl_vector_set(best,2,1.0);
  }

  for (int i=1;i<=10;i++){
    iter = 0;

    for (int zi=0;zi<znum;zi++){
      gsl_vector_set(x,0+zi*3,zmean[zi]);
      gsl_vector_set(x,1+zi*3,zsigma[zi]);
      gsl_vector_set(ss,0+zi*3,zmean[zi]/4.0);
      gsl_vector_set(ss,1+zi*3,zsigma[zi]/4.0);
      if(zi > 0){
	gsl_vector_set(x,2,1.0);
	gsl_vector_set(ss,2,0.1);
      }
    }
    
    s = gsl_multimin_fminimizer_alloc(T,vecsize);
    gsl_multimin_fminimizer_set(s,&minex_func,x,ss);
    
    printf("%5s \t %5s \t %5s","Trial","Mean","Sigma");
    if(vecsize == 5)
      printf(" \t %5s \t %5s \t %5s","Ratio","Mean","Sigma");
    printf(" \t %5s\n","Chi-Square");
    printf("%s\n","_______________________________________________________________");

    do{
      iter++;
      status = gsl_multimin_fminimizer_iterate(s);
      
      if(status)
	break;
      
      size = gsl_multimin_fminimizer_size(s);
      status = gsl_multimin_test_size(size,1e-2);
      
      if(status == GSL_SUCCESS){
	printf("%s\n","-----------------------------------------------------------");
	printf("\nFit %u Converged At Iteration %u\n",i,int(iter));
	printf("\nChi-Square: %7.3f\n\n",chi_min);
	printf("%s\n","Redshift Distribution 1: ");
	printf("   Mean:  %5.3f\n",gsl_vector_get(best,0));
	printf("   Sigma: %5.3f\n",abs(gsl_vector_get(best,1)));
	if(vecsize == 5){
	  printf("%s\n","Redshift Distribution 2: ");
	  printf("   Size Ratio: %5.3f\n",abs(gsl_vector_get(best,2)));
	  printf("   Mean:       %5.3f\n",gsl_vector_get(best,3));
	  printf("   Sigma:      %5.3f\n",abs(gsl_vector_get(best,4)));
	}

	if((s->fval < chi_min) or (i == 1)){
	  for (int bi=0;bi<vecsize;bi++)
	    gsl_vector_set(best,bi,gsl_vector_get(s->x,bi));
	  chi_min = s->fval;
	  if(tester.save(outfile))
	    printf("%s\n","Result Saved");
	  else
	    printf("%s\n","Save Unsuccessful");
	}
      }
      else{
	printf("%3u \t %5.3f \t %5.3f \t ",int(iter),gsl_vector_get(s->x,0),abs(gsl_vector_get(s->x,1)));
	if(vecsize == 5)
	  printf("%5.3f \t %5.3f \t %5.3f \t ",abs(gsl_vector_get(s->x,2)),gsl_vector_get(s->x,3),abs(gsl_vector_get(s->x,4)));
	printf("%5.1f \n",s->fval);
      }
    }while( status == GSL_CONTINUE and iter < 100);
   
    printf("%s\n","_______________________________________________________________");

    gsl_multimin_fminimizer_free(s);
  }
  
  gsl_vector_free(x);
  gsl_vector_free(ss);

  printf("%s\n","_______________________________________________________________");
  printf("\nBest Chi-Square: %7.3f\n\n",chi_min);
  printf("%s\n","Redshift Distribution 1: ");
  printf("   Mean:  %5.3f\n",gsl_vector_get(best,0));
  printf("   Sigma: %5.3f\n",abs(gsl_vector_get(best,1)));
  if(vecsize == 5){
    printf("%s\n","Redshift Distribution 2: ");
    printf("   Size Ratio: %5.3f\n",abs(gsl_vector_get(best,2)));
    printf("   Mean:       %5.3f\n",gsl_vector_get(best,3));
    printf("   Sigma:      %5.3f\n",abs(gsl_vector_get(best,4)));
  }
  printf("%s\n","_______________________________________________________________");

  gsl_vector_free(best);
  delete[] ps.p;

  return 0;
}

double simulate(const gsl_vector *v,void *params){
  static fixed_params *p;
  p = (fixed_params *)params; 
  return p->sim->simulate(v,params);
}
