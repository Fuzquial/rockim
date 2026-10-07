#include "rockim/PotentialContact.hpp"
#include <cstdio>
using namespace rockim; using pot::V2;
// generalized coords: 12 = 6 nodes x 2
static void forces(const double* q, double p, double* f, double* area){
  V2 a[3], b[3]; for(int k=0;k<3;++k){a[k]=V2(q[2*k],q[2*k+1]); b[k]=V2(q[6+2*k],q[6+2*k+1]);}
  pot::PairForce R; for(int i=0;i<12;++i) f[i]=0; *area=0;
  if(!pot::pairForce(a,b,p,R)) return;
  for(int k=0;k<3;++k){f[2*k]=R.fA[k].x();f[2*k+1]=R.fA[k].y();f[6+2*k]=R.fB[k].x();f[6+2*k+1]=R.fB[k].y();}
  *area=R.area;
}
// loop work along closed path q(s)=q0+ r*(cos s * u + sin s * w), midpoint rule with N steps
double loopWork(const double* q0,const double* u,const double* w,double r,double p,int N){
  double W=0, q[12], f[12], ar;
  for(int n=0;n<N;++n){ double s=2*M_PI*(n+0.5)/N, ds=2*M_PI/N;
    for(int i=0;i<12;++i) q[i]=q0[i]+r*(cos(s)*u[i]+sin(s)*w[i]);
    forces(q,p,f,&ar);
    for(int i=0;i<12;++i) W+=f[i]*r*(-sin(s)*u[i]+cos(s)*w[i])*ds; }
  return W;
}
int main(){
  double p=1.0;
  // A: triangle (0,0),(1,0),(0,1) CCW ; B: triangle shifted, overlapping on hypotenuse region
  double q0[12]={0,0, 1,0, 0,1,   0.9,0.9-0.05, 0.9-1,0.9-0.05+0.0, 0.9,0.9-0.05-1};
  // make B CCW: (0.9,0.85),( -0.1,0.85)?? check orientation below
  V2 b0(q0[6],q0[7]),b1(q0[8],q0[9]),b2(q0[10],q0[11]);
  double o=pot::cross2(b1-b0,b2-b0); printf("orientB %g\n",o);
  if(o<0){ std::swap(q0[8],q0[10]); std::swap(q0[9],q0[11]); }
  double f[12],ar; forces(q0,p,f,&ar); printf("area %g  F_A=(%g,%g)\n",ar,f[0]+f[2]+f[4],f[1]+f[3]+f[5]);
  double sumx=0,sumy=0,mom=0; for(int k=0;k<6;++k){sumx+=f[2*k];sumy+=f[2*k+1]; mom+=q0[2*k]*f[2*k+1]-q0[2*k+1]*f[2*k];}
  printf("sum F (%g,%g) moment %g\n",sumx,sumy,mom);
  double r=0.02; int N=20000;
  // 1) rigid translation of A in x,y
  double u[12]={0},w[12]={0}; for(int k=0;k<3;++k){u[2*k]=1;w[2*k+1]=1;}
  printf("loop rigid translation A        : W=%+.3e (F scale %g, r %g)\n",loopWork(q0,u,w,r,p,N),std::hypot(f[0]+f[2]+f[4],f[1]+f[3]+f[5]),r);
  // 2) translation x + rotation of A about its centroid (linearized rotation)
  double cx=1.0/3,cy=1.0/3; double u2[12]={0},w2[12]={0};
  for(int k=0;k<3;++k){u2[2*k]=1; w2[2*k]=-(q0[2*k+1]-cy); w2[2*k+1]=(q0[2*k]-cx);}
  printf("loop translation x + rotation A : W=%+.3e\n",loopWork(q0,u2,w2,r,p,N));
  // 3) single node loops: each node of A and B in x/y
  for(int nd=0;nd<6;++nd){ double u3[12]={0},w3[12]={0}; u3[2*nd]=1; w3[2*nd+1]=1;
    printf("loop node %d (x,y)               : W=%+.3e\n",nd,loopWork(q0,u3,w3,r,p,N)); }
  // 4) mixed: node0 x with node 4 y
  for(int a=0;a<12;++a) for(int b=a+1;b<12;++b){ double u4[12]={0},w4[12]={0}; u4[a]=1; w4[b]=1;
     double W=loopWork(q0,u4,w4,r,p,N); if(fabs(W)>1e-9) printf("  dof %2d,%2d : W=%+.3e  curl=%+.3e\n",a,b,W,W/(M_PI*r*r)); }
  // reference : energy-like scale p*area
  printf("p*area = %g\n",p*ar);
}
