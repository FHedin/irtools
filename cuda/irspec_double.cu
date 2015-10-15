#ifdef USE_DOUBLE
#define real_t double
#define fftComplex_t cufftDoubleComplex
#define complex_t cuDoubleComplex
#else
#define real_t double
#define fftComplex_t cufftDoubleComplex
#define complex_t cuDoubleComplex
#endif

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <cuda.h>
#include <cuda_runtime.h>
#include <cufft.h>
#include <cuComplex.h>

#define nThrdsX 16
#define nThrdsY 16
#define nThrds 256

#define clight 299792458
#define PI 3.14159265358979323846
#define cmtops (2*PI*clight*1.e-10)

__device__ __forceinline__ cuDoubleComplex my_cexpc (cuDoubleComplex z)
{

  cuDoubleComplex res;

  double t = exp(z.x);

  res.x=cos(z.y);
  res.y=sin(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__device__ __forceinline__ cuDoubleComplex my_cexpf (double z)
{

  cuDoubleComplex res;

  res.x=exp(z);
  res.y=0.0;

  return res;

}

__device__ __forceinline__ cuDoubleComplex my_cexpi (double z)
{

  cuDoubleComplex res;

  res.x=cos(z);
  res.y=sin(z);

  return res;

}

__device__ double g(double t, double *d_param)
{
  
  double f1,f2,f3,f4;
  
  f1=d_param[3]*exp(-t/d_param[6]);
  f2=d_param[4]*exp(-t/d_param[7]);
  f3=d_param[5]*exp(-t/d_param[8]);
  f4=(d_param[0]+d_param[1]+d_param[2])*t-(d_param[3]+d_param[4]+d_param[5]);
  
  return ( f1 + f2 + f3 + f4 );
  
}

__device__ cuDoubleComplex Rr(double t1 ,double t2, double t3, double *d_param)
{
  cuDoubleComplex f1;
  cuDoubleComplex f2;
  
  f1 = my_cexpi( -( t3 - t1 ) * d_param[9] * cmtops );
  f1 = cuCsub(f1 , cuCmul( my_cexpi( -cmtops * ( ( d_param[9] - d_param[10] ) * t3 - ( d_param[9] * t1 ) ) ) , my_cexpf( -(d_param[12]*t3)/(2.*d_param[11]) ) ) );
  f1 = cuCmul(f1 , my_cexpf( -( t1 + t3 + (2.*t2) ) / ( 2.*d_param[11] ) ) );
  f2 = my_cexpf(-g(t1,d_param) + g(t2,d_param) - g(t3,d_param) - g(t1+t2,d_param) - g(t2+t3,d_param) + g(t1+t2+t3,d_param)  );
  
  return cuCmul( f1 , f2 );
}

__device__ cuDoubleComplex Rnr(double t1 ,double t2, double t3, double *d_param)
{
  cuDoubleComplex f1;
  cuDoubleComplex f2;
  
  f1 = my_cexpi( -( t3 + t1 ) * d_param[9] * cmtops );
  f1 = cuCsub(f1 , cuCmul( my_cexpi( -cmtops * ( ( d_param[9] - d_param[10] ) * t3 + ( d_param[9] * t1 ) ) ) , my_cexpf( -(d_param[12]*t3)/(2.*d_param[11]) ) ) );
  f1 = cuCmul(f1 , my_cexpf( -( t1 + t3 + (2.*t2) ) / ( 2.*d_param[11] ) ) );
  f2 = my_cexpf( -g(t1,d_param) - g(t2,d_param) - g(t3,d_param) + g(t1+t2,d_param) + g(t2+t3,d_param) - g(t1+t2+t3,d_param) );
  
  return cuCmul( f1 , f2 );
}

__global__ void kernelHalfFirstX(cufftDoubleComplex *d_r)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  cuDoubleComplex half;
  
  half.x=0.5;
  half.y=0.0;
  
  d_r[i]=cuCmul(d_r[i],half);
  
}

__global__ void kernelHalfFirstY(cufftDoubleComplex *d_r, int n)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  int j = i*n;
  
  cuDoubleComplex half;
  
  half.x=0.5;
  half.y=0.0;
  
  if (j>0)
    d_r[j]=cuCmul(d_r[j],half); 
  
}

__global__ void kernelPermutation(cufftDoubleComplex *d_r,cufftDoubleComplex *d_ftr, int n)
{
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  int k = i * n + j;
  int l = i * n + j + n/2 ;
  int m = ( i + n/2 ) * n + j;
  int p = ( i + n/2 ) * n + j + n/2;
  
  d_r[m]=d_ftr[l];
  d_r[p]=d_ftr[k];
  d_r[l]=d_ftr[m];
  d_r[k]=d_ftr[p];
}

__global__ void kernelRephasing(cufftDoubleComplex *d_r,double dt,double t2,int n,double *d_param)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  int k = i * n + j;
  
  double t1 = i * dt;
  double t3 = j * dt;
  
  d_r[k]=Rr(t1,t2,t3,d_param);
  
}

__global__ void kernelNonRephasing(cufftDoubleComplex *d_r, double dt, double t2,int n,double *d_param)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  int k= i * n + j;
  
  double t1 = i * dt;
  double t3 = j * dt;
  
  d_r[k]=Rnr(t1,t2,t3,d_param);
  
}

int main(int argc, char* argv[])
{
  
  FILE* input=fopen(argv[1],"r");
  FILE* out=fopen("spec.dat","w");
  
  int i,ii,j,l,test;
  int iw1,iw3,iFirst,iLast,jFirst,jLast;
  int ndata,nw1,nw3,nave;
  
  double t2,dt;
  double w1min,w1max,w3min,w3max;
  double w1,w3;
  double **res;
  
  double c1,tau1,c2,tau2,c3,tau3,wm,delta,tLife,alpha;
  
  fscanf(input,"%lf %lf %lf %lf",&w1min,&w1max,&w3min,&w3max);
  fscanf(input,"%lf %lf %lf",&delta,&tLife,&alpha);
  fscanf(input,"%d %lf",&ndata,&dt);
  fscanf(input,"%lf",&t2);
  fscanf(input,"%d",&nave);
  
  ndata=(int)(ndata/nThrds)+1;
  ndata=nThrds*2*(int)(ndata/2);
 
  printf("%d %d\n",ndata,nThrds);
 
  nw1=0;
  test=1;
  for(i=ndata/2-1;i<ndata-1;i++)
  {
    
    w1=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+i);
    
    if(w1>w1max)
    {
      iLast=i-1;
      break;
    }
    
    if( w1>=w1min )
    {
      nw1++;
      
      if(test)
      {
	test=0;
	iFirst=i;
      }
    }
    
  }
  
  nw3=0;
  test=1;
  for(j=ndata/2-1;j<ndata-1;j++)
  {
    
    w3=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+j);
    
    if(w3>w3max)
    {
      jLast=j-1;
      break;
    }
    
    if( w3>=w3min )
    {
      nw3++;
      
      if(test)
      {
	test=0;
	jFirst=j;
      }
    }
    
  }

  printf("%d %d %d %d %d %d\n",iFirst,iLast,jFirst,jLast,nw1,nw3);
  
  res=(double**)malloc(nw1*sizeof(double*));
  
  for(i=0;i<nw1;i++)
  {
    res[i]=(double*)malloc(nw3*sizeof(double));
    for(j=0;j<nw3;j++)
      res[i][j]=0.0;
  }
  
  double *param,*d_param;
  
  param=(double*)malloc(13*sizeof(double));
  cudaMalloc((void**) &d_param,13*sizeof(double));
  
  param[10]=delta;
  param[11]=tLife;
  param[12]=alpha;
  
  cufftDoubleComplex *d_rlist,*d_ftrlist,*rlist;
  cufftHandle d_pr;
  
  rlist=(cufftDoubleComplex*)malloc(ndata*ndata*sizeof(cufftDoubleComplex));
  
  cudaMalloc((void**) &d_rlist,ndata*ndata*sizeof(cufftDoubleComplex));
  cudaMalloc((void**) &d_ftrlist,ndata*ndata*sizeof(cufftDoubleComplex));
  
  cufftPlan2d(&d_pr,ndata,ndata,CUFFT_Z2Z);
  
  dim3 threadsPerBlock(nThrdsX, nThrdsY);
  dim3 numBlocks(ndata / threadsPerBlock.x, ndata / threadsPerBlock.y);
  dim3 numHalfBlocks( (ndata/2) / threadsPerBlock.x, (ndata/2) / threadsPerBlock.y );
  
  for(ii=0;ii<nave;ii++)
  {
    
    fscanf(input,"%lf",&wm);
    fscanf(input,"%lf",&c1);
    fscanf(input,"%lf",&c2);
    fscanf(input,"%lf",&c3);
    fscanf(input,"%lf",&tau1);
    fscanf(input,"%lf",&tau2);
    fscanf(input,"%lf",&tau3);
    
    printf("%lf %lf %lf %lf %lf %lf %lf\n",wm,c1,c2,c3,tau1,tau2,tau3);

    param[0]=c1*tau1;
    param[1]=c2*tau2;
    param[2]=c3*tau3;
    
    param[3]=param[0]*tau1;
    param[4]=param[1]*tau2;
    param[5]=param[2]*tau3;
    
    param[6]=tau1;
    param[7]=tau2;
    param[8]=tau3;
    param[9]=wm;
    
    cudaMemcpy(d_param,param,13*sizeof(double),cudaMemcpyHostToDevice);
    
    // Rephasing component
    
    kernelRephasing<<<numBlocks,threadsPerBlock>>>(d_rlist,dt,t2,ndata,d_param);
    kernelHalfFirstX<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist);
    kernelHalfFirstY<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist,ndata);
    
    cufftExecZ2Z(d_pr,d_rlist,d_ftrlist,CUFFT_INVERSE);
    
    kernelPermutation<<<numHalfBlocks,threadsPerBlock>>>(d_rlist,d_ftrlist,ndata);
    
    cudaMemcpy(rlist,d_rlist,ndata*ndata*sizeof(cufftDoubleComplex),cudaMemcpyDeviceToHost);
    
    iw1=0;
    for(i=iFirst;i<=iLast;i++)
    {
      iw3=0;
      for(j=jFirst;j<=jLast;j++)
      {
	
	l=ndata*(ndata-i-1)+j+1;
	res[iw1][iw3]+=cuCreal(rlist[l]);
	
	iw3++;
      }
      iw1++;
    }
    
    // Non-rephasing component
    
    kernelNonRephasing<<<numBlocks,threadsPerBlock>>>(d_rlist,dt,t2,ndata,d_param);
    kernelHalfFirstX<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist);
    kernelHalfFirstY<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist,ndata);
    
    cufftExecZ2Z(d_pr,d_rlist,d_ftrlist,CUFFT_INVERSE);
    
    kernelPermutation<<<numHalfBlocks,threadsPerBlock>>>(d_rlist,d_ftrlist,ndata);
    
    cudaMemcpy(rlist,d_rlist,ndata*ndata*sizeof(cufftDoubleComplex),cudaMemcpyDeviceToHost);
    
    iw1=0;
    for(i=iFirst;i<=iLast;i++)
    {
      iw3=0;
      for(j=jFirst;j<=jLast;j++)
      {
	
	l=ndata*(i+1)+j+1;
	res[iw1][iw3]+=cuCreal(rlist[l]);
	
	iw3++;
      }
      iw1++;
    }
  }
  
  cufftDestroy(d_pr);
  
  cudaFree(d_rlist);
  cudaFree(d_ftrlist);
  
  iw1=0;
  for(i=iFirst;i<=iLast;i++)
  {
    
    w1=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+i);
    
    iw3=0;
    for(j=jFirst;j<=jLast;j++)
    {
      
      w3=1.e10/((double)ndata*clight*dt)*(double)(1-ndata/2+j);
      
      fprintf(out,"%lf\t%lf\t%le\n",w1,w3,res[iw1][iw3]/(double)ndata);
      iw3++;
    }
    
    fprintf(out,"\n");
    iw1++;
  }
  
  for(i=0;i<nw1;i++)
    free(res[i]);
  
  free(res);
  
  fclose(out);
  fclose(input);
  
  return 0;
  
}
