#include "potx.hpp"
#include <cstdio>
using namespace rockim; using pot::V2;
bool corr=false;
static double energy(const double*q,double p){ V2 a[3],b[3]; for(int k=0;k<3;++k){a[k]=V2(q[2*k],q[2*k+1]);b[k]=V2(q[6+2*k],q[6+2*k+1]);}
  double IA,IB,ar; enq::integrals(a,b,IA,IB,ar); return p*(IA+IB);}
static void forces(const double* q, double p, double* f){
  V2 a[3], b[3]; for(int k=0;k<3;++k){a[k]=V2(q[2*k],q[2*k+1]); b[k]=V2(q[6+2*k],q[6+2*k+1]);}
  pot::PairForce R; for(int i=0;i<12;++i) f[i]=0;
  if(!pot::pairForce(a,b,p,R)) return; if(corr) enq::correct(a,b,p,R);
  for(int k=0;k<3;++k){f[2*k]=R.fA[k].x();f[2*k+1]=R.fA[k].y();f[6+2*k]=R.fB[k].x();f[6+2*k+1]=R.fB[k].y();}}
double loopWork(const double* q0,const double* u,const double* w,double r,double p,int N){
  double W=0,q[12],f[12]; for(int n=0;n<N;++n){double s=2*M_PI*(n+0.5)/N,ds=2*M_PI/N;
    for(int i=0;i<12;++i) q[i]=q0[i]+r*(cos(s)*u[i]+sin(s)*w[i]); forces(q,p,f);
    for(int i=0;i<12;++i) W+=f[i]*r*(-sin(s)*u[i]+cos(s)*w[i])*ds;} return W;}
int main(){
  double p=1.0; double q0[12]={0,0, 1,0, 0,1,  0.9,0.85, 0.9,-0.15, -0.1,0.85};
  // orientation of B
  V2 b0(q0[6],q0[7]),b1(q0[8],q0[9]),b2(q0[10],q0[11]); if(pot::cross2(b1-b0,b2-b0)<0){std::swap(q0[8],q0[10]);std::swap(q0[9],q0[11]);}
  // gradient check : -dE/dq vs Munjiza and corrected
  double fM[12],fC[12]; corr=false; forces(q0,p,fM); corr=true; forces(q0,p,fC);
  printf(" dof   -dE/dq(FD)     Munjiza      corrige\n");
  double eM=0,eC=0,nn=0;
  for(int i=0;i<12;++i){ double qp[12],qm[12]; for(int j=0;j<12;++j){qp[j]=qm[j]=q0[j];} double h=1e-6; qp[i]+=h; qm[i]-=h;
    double g=-(energy(qp,p)-energy(qm,p))/(2*h); printf(" %2d  %+.6f  %+.6f  %+.6f\n",i,g,fM[i],fC[i]);
    eM=std::max(eM,fabs(fM[i]-g)); eC=std::max(eC,fabs(fC[i]-g)); nn=std::max(nn,fabs(g)); }
  printf("ecart max / |f|max : Munjiza %.3e   corrige %.3e\n",eM/nn,eC/nn);
  for(int c=0;c<2;++c){ corr=c; double wmax=0; int ia=0,ib=0;
    for(int a=0;a<12;++a) for(int b=a+1;b<12;++b){ double u[12]={0},w[12]={0}; u[a]=1; w[b]=1; double W=fabs(loopWork(q0,u,w,0.02,p,4000)); if(W>wmax){wmax=W;ia=a;ib=b;} }
    printf("%s : max |travail sur boucle| = %.3e (ddl %d,%d ; r=0.02)\n",c?"corrige":"Munjiza",wmax,ia,ib); }
}
