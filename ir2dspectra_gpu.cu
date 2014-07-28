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

__global__ void kernelHalfFirstX(cufftComplex* d_r);
__global__ void kernelHalfFirstY(cufftComplex* d_r, int n);
__global__ void kernelPermutation(cufftComplex* d_r,cufftComplex* d_ftr, int n);
__global__ void kernelRephasing(cufftComplex* d_r,float dt,float t2,int n,float* d_param);
__global__ void kernelNonRephasing(cufftComplex* d_r, float dt, float t2,int n,float* d_param);

__device__ float g(float t,float* d_param);
__device__ cuComplex Rr(float t1 ,float t2, float t3,float* d_param);
__device__ cuComplex Rnr(float t1 ,float t2, float t3,float* d_param);

__device__ __forceinline__ cuComplex my_cexpc (cuComplex z);
__device__ __forceinline__ cuComplex my_cexpf (float z);
__device__ __forceinline__ cuComplex my_cexpi (float z);

int main(int argc, char* argv[])
{
  
  FILE* input=fopen(argv[1],"r");
  FILE* out=fopen("spec.dat","w");
  
  int i,ii,j,l,test;
  int iw1,iw3,iFirst,iLast,jFirst,jLast;
  int ndata,nw1,nw3,nave;
  
  float t2,dt;
  float w1min,w1max,w3min,w3max;
  float w1,w3;
  float **res;
  
  float c1,tau1,c2,tau2,c3,tau3,wm,delta,tLife,alpha;
  
  fscanf(input,"%f %f %f %f",&w1min,&w1max,&w3min,&w3max);
  fscanf(input,"%f %f %f",&delta,&tLife,&alpha);
  fscanf(input,"%d %f",&ndata,&dt);
  fscanf(input,"%f",&t2);
  fscanf(input,"%d",&nave);
  
  ndata=(int)(ndata/nThrds)+1;
  ndata=nThrds*2*(int)(ndata/2);
  
  nw1=0;
  test=1;
  for(i=ndata/2-1;i<ndata-1;i++)
  {
    
    w1=1.e10/((float)ndata*clight*dt)*(float)(1-ndata/2+i);
    
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
    
    w3=1.e10/((float)ndata*clight*dt)*(float)(1-ndata/2+j);
    
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
  
  res=(float**)malloc(nw1*sizeof(float*));
  
  for(i=0;i<nw1;i++)
  {
    res[i]=(float*)malloc(nw3*sizeof(float));
    for(j=0;j<nw3;j++)
      res[i][j]=0.0f;
  }
  
  float *param,*d_param;
  
  param=(float*)malloc(13*sizeof(float));
  cudaMalloc((void**) &d_param,13*sizeof(float));
  
  param[10]=delta;
  param[11]=tLife;
  param[12]=alpha;
  
  cufftComplex *d_rlist,*d_ftrlist,*rlist;
  cufftHandle d_pr;
  
  rlist=(cufftComplex*)malloc(ndata*ndata*sizeof(cufftComplex));
  
  cudaMalloc((void**) &d_rlist,ndata*ndata*sizeof(cufftComplex));
  cudaMalloc((void**) &d_ftrlist,ndata*ndata*sizeof(cufftComplex));
  
  cufftPlan2d(&d_pr,ndata,ndata,CUFFT_C2C);
  
  dim3 threadsPerBlock(nThrdsX, nThrdsY);
  dim3 numBlocks(ndata / threadsPerBlock.x, ndata / threadsPerBlock.y);
  dim3 numHalfBlocks( (ndata/2) / threadsPerBlock.x, (ndata/2) / threadsPerBlock.y );
  
  for(ii=0;ii<nave;ii++)
  {
    
    fscanf(input,"%f",&wm);
    fscanf(input,"%f",&c1);
    fscanf(input,"%f",&c2);
    fscanf(input,"%f",&c3);
    fscanf(input,"%f",&tau1);
    fscanf(input,"%f",&tau2);
    fscanf(input,"%f",&tau3);
    
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
    
    cudaMemcpy(d_param,d_param,13*sizeof(float),cudaMemcpyHostToDevice);
    
    // Rephasing component
    
    kernelRephasing<<<numBlocks,threadsPerBlock>>>(d_rlist,dt,t2,ndata,d_param);
    kernelHalfFirstX<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist);
    kernelHalfFirstY<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist,ndata);
    
    cufftExecC2C(d_pr,d_rlist,d_ftrlist,CUFFT_INVERSE);
    
    kernelPermutation<<<numHalfBlocks,threadsPerBlock>>>(d_rlist,d_ftrlist,ndata);
    
    cudaMemcpy(rlist,d_rlist,ndata*ndata*sizeof(cufftComplex),cudaMemcpyDeviceToHost);
    
    iw1=0;
    for(i=iFirst;i<=iLast;i++)
    {
      iw3=0;
      for(j=jFirst;j<=jLast;j++)
      {
	
	l=ndata*(ndata-i-1)+j+1;
	res[iw1][iw3]=cuCrealf(rlist[l]);
	
	iw3++;
      }
      iw1++;
    }
    
    // Non-rephasing component
    
    kernelNonRephasing<<<numBlocks,threadsPerBlock>>>(d_rlist,dt,t2,ndata,d_param);
    kernelHalfFirstX<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist);
    kernelHalfFirstY<<<(ndata / threadsPerBlock.x),threadsPerBlock>>>(d_rlist,ndata);
    
    cufftExecC2C(d_pr,d_rlist,d_ftrlist,CUFFT_INVERSE);
    
    kernelPermutation<<<numHalfBlocks,threadsPerBlock>>>(d_rlist,d_ftrlist,ndata);
    
    cudaMemcpy(rlist,d_rlist,ndata*ndata*sizeof(cufftComplex),cudaMemcpyDeviceToHost);
    
    iw1=0;
    for(i=iFirst;i<=iLast;i++)
    {
      iw3=0;
      for(j=jFirst;j<=jLast;j++)
      {
	
	l=ndata*(i+1)+j+1;
	res[iw1][iw3]=cuCrealf(rlist[l]);
	
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
    
    w1=1.e10/((float)ndata*clight*dt)*(float)(1-ndata/2+i);
    
    iw3=0;
    for(j=jFirst;j<=jLast;j++)
    {
      
      w3=1.e10/((float)ndata*clight*dt)*(float)(1-ndata/2+j);
      
      fprintf(out,"%f\t%f\t%e\n",w1,w3,res[iw1][iw3]/(float)ndata);
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

__global__ void kernelHalfFirstX(cufftComplex* d_r)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  cuComplex half;
  
  half.x=0.5f;
  half.y=0.0f;
  
  d_r[i]=cuCmulf(d_r[i],half);
  
}

__global__ void kernelHalfFirstY(cufftComplex* d_r, int n)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  
  int j = i*n;
  
  cuComplex half;
  
  half.x=0.5f;
  half.y=0.0f;
  
  if (j>0)
    d_r[j]=cuCmulf(d_r[j],half); 
  
}

__global__ void kernelPermutation(cufftComplex* d_r,cufftComplex* d_ftr, int n)
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

__global__ void kernelRephasing(cufftComplex* d_r,float dt,float t2,int n,float* d_param)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  int k = i * n + j;
  
  float t1 = i * dt;
  float t3 = j * dt;
  
  d_r[k]=Rr(t1,t2,t3,d_param);
  
}

__global__ void kernelNonRephasing(cufftComplex* d_r, float dt, float t2,int n,float* d_param)
{
  
  int i = blockIdx.x * blockDim.x + threadIdx.x;
  int j = blockIdx.y * blockDim.y + threadIdx.y;
  
  int k= i * n + j;
  
  float t1 = i * dt;
  float t3 = j * dt;
  
  d_r[k]=Rnr(t1,t2,t3,d_param);
  
}

__device__ float g(float t, float* d_param)
{
  
  float f1,f2,f3,f4;
  
  f1=d_param[3]*expf(-t/d_param[6]);
  f2=d_param[4]*expf(-t/d_param[7]);
  f3=d_param[5]*expf(-t/d_param[8]);
  f4=(d_param[0]+d_param[1]+d_param[2])*t-(d_param[3]+d_param[4]+d_param[5]);
  
  return ( f1 + f2 + f3 + f4 );
  
}

__device__ cuComplex Rr(float t1 ,float t2, float t3, float* d_param)
{
  cuComplex f1;
  cuComplex f2;
  
  f1 = my_cexpi( -( t3 - t1 ) * d_param[9] * cmtops );
  f1 = cuCsubf(f1 , cuCmulf( my_cexpi( -cmtops * ( ( d_param[9] - d_param[10] ) * t3 - ( d_param[9] * t1 ) ) ) , my_cexpf( -(d_param[12]*t3)/(2.*d_param[11]) ) ) );
  f1 = cuCmulf(f1 , my_cexpf( -( t1 + t3 + (2.*t2) ) / ( 2.*d_param[11] ) ) );
  f2 = my_cexpf(-g(t1,d_param) + g(t2,d_param) - g(t3,d_param) - g(t1+t2,d_param) - g(t2+t3,d_param) + g(t1+t2+t3,d_param)  );
  
  return cuCmulf( f1 , f2 );
}

__device__ cuComplex Rnr(float t1 ,float t2, float t3, float* d_param)
{
  cuComplex f1;
  cuComplex f2;
  
  f1 = my_cexpi( -( t3 + t1 ) * d_param[9] * cmtops );
  f1 = cuCsubf(f1 , cuCmulf( my_cexpi( -cmtops * ( ( d_param[9] - d_param[10] ) * t3 + ( d_param[9] * t1 ) ) ) , my_cexpf( -(d_param[12]*t3)/(2.*d_param[11]) ) ) );
  f1 = cuCmulf(f1 , my_cexpf( -( t1 + t3 + (2.*t2) ) / ( 2.*d_param[11] ) ) );
  f2 = my_cexpf( -g(t1,d_param) - g(t2,d_param) - g(t3,d_param) + g(t1+t2,d_param) + g(t2+t3,d_param) - g(t1+t2+t3,d_param) );
  
  return cuCmulf( f1 , f2 );
}

__device__ __forceinline__ cuComplex my_cexpc (cuComplex z)
{

  cuComplex res;

  float t = expf(z.x);

  res.x=cosf(z.y);
  res.y=sinf(z.y);

  res.x *= t;
  res.y *= t;

  return res;

}

__device__ __forceinline__ cuComplex my_cexpf (float z)
{

  cuComplex res;

  res.x=expf(z);
  res.y=0.f;

  return res;

}

__device__ __forceinline__ cuComplex my_cexpi (float z)
{

  cuComplex res;

  res.x=cosf(z);
  res.y=sinf(z);

  return res;

}
