#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <fftw3.h>

#define clight 299792458
#define PI 3.14159265358979323846

double c1,tau1,c2,tau2,c3,tau3,w,wm,delta;
double cmtops=2*PI*clight*1.e-10;

double g(double t);
double complex Rr(double t1 ,double t2, double t3);
double complex Rnr(double t1 ,double t2, double t3);

int main(int argc, char* argv[])
{
  FILE* input=fopen(argv[1],"r");
  FILE* out=fopen("spec.dat","w");
  int i,j,k,kk,l,ll,ndata;
  double t1,t2,t3,dt;
  double w1min,w1max,w3min,w3max;
  
  fftw_complex *rlist,*ftrlist,*rnlist,*ftrnlist;
  fftw_plan pr,prn;
  
  fscanf(input,"%lf",&w);
  fscanf(input,"%lf",&c1);
  fscanf(input,"%lf",&c2);
  fscanf(input,"%lf",&c3);
  fscanf(input,"%lf",&tau1);
  fscanf(input,"%lf",&tau2);
  fscanf(input,"%lf",&tau3);
  fscanf(input,"%lf",&delta);
  fscanf(input,"%lf",&wm);
  fscanf(input,"%d",&ndata);
  fscanf(input,"%lf",&dt);
  fscanf(input,"%lf %lf %lf %lf",&w1min,&w1max,&w3min,&w3max);
  fscanf(input,"%lf",&t2);
  
  rlist=(fftw_complex*) fftw_malloc(ndata*ndata*sizeof(*rlist));
  ftrlist=(fftw_complex*) fftw_malloc(ndata*ndata*sizeof(*ftrlist));
  rnlist=(fftw_complex*) fftw_malloc(ndata*ndata*sizeof(*rnlist));
  ftrnlist=(fftw_complex*) fftw_malloc(ndata*ndata*sizeof(*ftrnlist));
  
  pr=fftw_plan_dft_2d(ndata, ndata, rlist,ftrlist , FFTW_BACKWARD, FFTW_ESTIMATE);
  prn=fftw_plan_dft_2d(ndata, ndata, rnlist,ftrnlist , FFTW_BACKWARD, FFTW_ESTIMATE);
    
  k=0;
  for(i=0;i<ndata;i++)
  {
    
    t1=(double)i*dt;
    
    for(j=0;j<ndata;j++)
    {
      
      t3=(double)j*dt;
      
      rlist[k]=Rr(t1,t2,t3);
      rnlist[k]=Rnr(t1,t2,t3);
      
      k++;
    }
    
  }
  
  for(i=0;i<ndata;i++)
  {
      
    rlist[i]=0.5*rlist[i];
    rnlist[i]=0.5*rnlist[i];
    
  }
  
  for(i=1;i<ndata;i++)
  {
    k=i*ndata;
    rlist[k]=0.5*rlist[k];
    rnlist[k]=0.5*rnlist[k];
    
  }
  
  fftw_execute(pr);
  fftw_execute(prn);
  
  for(i=0;i<ndata/2;i++)
  {
    
    for(j=0;j<ndata/2;j++)
    {
        
      k=i*ndata+j;
      kk=i*ndata+j+ndata/2;
      l=(i+ndata/2)*ndata+j;
      ll=(i+ndata/2)*ndata+j+ndata/2;
      
      rlist[l]=ftrlist[kk];
      rlist[ll]=ftrlist[k];
      rlist[k]=ftrlist[ll];
      rlist[kk]=ftrlist[l];
      
      rnlist[l]=ftrnlist[kk];
      rnlist[ll]=ftrnlist[k];
      rnlist[k]=ftrnlist[ll];
      rnlist[kk]=ftrnlist[l];
      
      
    }
    
  }
  
  fftw_free(ftrlist);
  fftw_free(ftrnlist);
  
  ftrlist=(fftw_complex*) fftw_malloc((ndata-1)*(ndata-1)*sizeof(*ftrlist));
  ftrnlist=(fftw_complex*) fftw_malloc((ndata-1)*(ndata-1)*sizeof(*ftrnlist));
  
  k=0;
  for(i=1;i<ndata;i++)
  {
    
    for(j=1;j<ndata;j++)
    {
      
      l=ndata*i+j;
      ll=ndata*(ndata-i)+j;
      
      ftrlist[k]=rlist[ll];
      ftrnlist[k]=rnlist[l];
      
      k++;
    }
    
  }
  
  double w1,w3,res;
  
  k=0;
  for(i=ndata/2-1;i<ndata-1;i++)
  {
    
    w1=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+i);
    
    if(w1>w1max)
      break;
    
    if( w1>=w1min )
    {
      
      for(j=ndata/2-1;j<ndata-1;j++)
      {
	
	k=i*(ndata-1)+j;
	w3=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+j);
	
	if(w3>w3max)
	  break;
	
	if( w3>=w3min )
	{
	  
	  res=creal(ftrlist[k]+ftrnlist[k])/(double)ndata;
	
	  fprintf(out,"%lf\t%lf\t%le\n",w1,w3,res);
	  
	}
	
      }
      
      fprintf(out,"\n");
      
    }
    
  }
  
  fftw_destroy_plan(pr);
  fftw_destroy_plan(prn);
  
  fftw_free(rlist);
  fftw_free(ftrlist);
  fftw_free(rnlist);
  fftw_free(ftrnlist);
  
  fclose(out);
  fclose(input);
  
  return 0;
  
}

double g(double t)
{
  
  double a,b,f1,f2,f3;
  
  a=w*tau1*tau1/(w*w*tau1*tau1+1.);
  b=tau1/(w*w*tau1*tau1+1.);
  
  f1=c1*exp(-t/tau1)*((b*b-a*a)*cos(w*t)-2.*a*b*sin(w*t));
  f2=c2*tau2*tau2*exp(-t/tau2)+c3*tau3*tau3*exp(-t/tau3);
  f3=(c1*b+c2*tau2+c3*tau3)*t-c1*(b*b-a*a)-(c2*tau2*tau2+c3*tau3*tau3);
  
  return ( f1 + f2 + f3 );
  
}

double complex Rr(double t1 ,double t2, double t3)
{
  double complex f1;
  double f2;
  
  f1 = cexp( -I * ( t3 - t1 ) * wm * cmtops ) - cexp( -I * cmtops * ( ( wm - delta ) * t3 - ( wm * t1 ) ) );
  f2 = exp( -g(t1) + g(t2) - g(t3) - g(t1+t2) - g(t2+t3) + g(t1+t2+t3) );
  
  return ( f1 * f2 );
}

double complex Rnr(double t1 ,double t2, double t3)
{
  double complex f1;
  double f2;
  
  f1 = cexp( -I * ( t3 + t1 ) * wm * cmtops ) - cexp( -I * cmtops * ( ( wm - delta ) * t3 + ( wm * t1 ) ) );
  f2 = exp( -g(t1) - g(t2) - g(t3) + g(t1+t2) + g(t2+t3) - g(t1+t2+t3) );
  
  return ( f1 * f2 );
}
