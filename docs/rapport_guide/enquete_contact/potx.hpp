// Outils d'enquete : energie exacte E = p int_S (phiA+phiB) dA, integrales I_A, I_B,
// et force corrigee (gradient exact) = Munjiza - p I_A grad(lambda_a) (A), idem B.
#pragma once
#include "rockim/PotentialContact.hpp"
namespace enq {
using rockim::pot::V2; using rockim::pot::cross2;
// clip polygon P (n) by half-plane a.x + b.y + c >= 0
inline int clipHP(const V2* P,int n,double a,double b,double c,V2* out){
  int m=0; for(int i=0;i<n;++i){const V2&p=P[i];const V2&q=P[(i+1)%n];
    double dp=a*p.x()+b*p.y()+c, dq=a*q.x()+b*q.y()+c;
    if(dp>=0){out[m++]=p; if(dq<0) out[m++]=p+(q-p)*(dp/(dp-dq));}
    else if(dq>=0) out[m++]=p+(q-p)*(dp/(dp-dq));
    if(m>14) break;}
  return m;}
// lambda_i(X) affine : coefficients g_i (gradient) and c_i so that l_i = g_i.X + c_i
struct Aff { V2 g[3]; double c[3]; V2 gradL[3];
  void set(const V2 P[3]){ double den=cross2(P[1]-P[0],P[2]-P[0]);
    for(int i=0;i<3;++i){ const V2&Pj=P[(i+1)%3]; const V2&Pk=P[(i+2)%3]; V2 e=Pk-Pj;
      // l_i = cross(Pk-Pj, X-Pj)/den
      g[i]=V2(-e.y(),e.x())/den; c[i]=-g[i].dot(Pj); gradL[i]=g[i]; } }
  double l(int i,const V2&X)const{return g[i].dot(X)+c[i];}
};
inline double polyArea(const V2*S,int n,V2&cen){double A2=0;cen.setZero();
  for(int i=0;i<n;++i){const V2&P=S[i];const V2&Q=S[(i+1)%n];double w=cross2(P,Q);A2+=w;cen+=(P+Q)*w;}
  if(std::fabs(A2)>1e-300) cen/=3*A2; return 0.5*A2;}
// integrals over S=A∩B of phiA and phiB
inline bool integrals(const V2 pa[3],const V2 pb[3],double&IA,double&IB,double&area){
  IA=IB=area=0; V2 S[8]; int n=rockim::pot::clipTriTri(pa,pb,S); if(n<3) return false;
  Aff fa,fb; fa.set(pa); fb.set(pb);
  for(int i=0;i<3;++i) for(int j=0;j<3;++j){
    V2 P1[16],P2[16]; int m=n; for(int k=0;k<n;++k)P1[k]=S[k];
    // region argmin_A = i : l_i - l_k <= 0 for k != i
    for(int k=0;k<3;++k){ if(k==i) continue; V2 g=fa.g[k]-fa.g[i]; double c=fa.c[k]-fa.c[i];
      m=clipHP(P1,m,g.x(),g.y(),c,P2); for(int t=0;t<m;++t)P1[t]=P2[t]; if(m<3)break;}
    if(m<3) continue;
    for(int k=0;k<3;++k){ if(k==j) continue; V2 g=fb.g[k]-fb.g[j]; double c=fb.c[k]-fb.c[j];
      m=clipHP(P1,m,g.x(),g.y(),c,P2); for(int t=0;t<m;++t)P1[t]=P2[t]; if(m<3)break;}
    if(m<3) continue;
    V2 cen; double a=polyArea(P1,m,cen); if(a<=0) continue;
    IA+=a*3*fa.l(i,cen); IB+=a*3*fb.l(j,cen); area+=a; }
  return true;}
// corrected (exact-gradient) nodal forces, in place on a PairForce from pairForce(.,.,p,R)
inline void correct(const V2 pa[3],const V2 pb[3],double p,rockim::pot::PairForce&R,double*IAo=nullptr,double*IBo=nullptr){
  double IA,IB,ar; integrals(pa,pb,IA,IB,ar); Aff fa,fb; fa.set(pa); fb.set(pb);
  for(int k=0;k<3;++k){ R.fA[k]-=p*IA*fa.gradL[k]; R.fB[k]-=p*IB*fb.gradL[k]; }
  if(IAo)*IAo=IA; if(IBo)*IBo=IB; }
}
