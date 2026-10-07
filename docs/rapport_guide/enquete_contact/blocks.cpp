// Banc autonome : deux blocs elastiques DEFORMABLES (CST, deformation plane) en contact
// par le potentiel de Munjiza (PotentialContact.hpp du depot, inchange), integres par le
// MEME schema que FdemSolver (forces a x_n, v += dt f/m, x += dt v ; compteur f.v_old dt).
// Energie mecanique exacte : KE + U_el + U_c, U_c = p int_S (phiA+phiB) dA (EnqPotx).
// usage: blocks mode(0 rebond,1 presse) corr(0 Munjiza,1 gradient exact) pfac dtfac n
#include "potx.hpp"
#include <cstdio>
#include <vector>
#include <random>
using namespace rockim; using pot::V2;
struct Tri{int n[3]; double A0; Eigen::Matrix<double,3,6> B; int body;};
int main(int argc,char**argv){
  int mode=atoi(argv[1]); int corr=atoi(argv[2]); double pfac=atof(argv[3]); double dtfac=atof(argv[4]); int N=argc>5?atoi(argv[5]):4;
  double E=1e10,nu=0.25,rho=2650,L=1.0; double p=pfac*E; // p = potPenaltyFactor * E * thk (thk=1)
  std::vector<V2> X; std::vector<double> m; std::vector<Tri> T; std::mt19937 rng(7); std::uniform_real_distribution<double> U(-0.15,0.15);
  auto block=[&](double x0,double y0,int body){ int base=X.size();
    for(int j=0;j<=N;++j) for(int i=0;i<=N;++i){ double px=x0+L*i/N, py=y0+L*j/N; bool edge=(i==0||j==0||i==N||j==N);
      if(!edge){px+=U(rng)*L/N; py+=U(rng)*L/N;} X.push_back(V2(px,py)); m.push_back(0);}
    for(int j=0;j<N;++j) for(int i=0;i<N;++i){ int a=base+j*(N+1)+i,b=a+1,c=a+N+1,d=c+1;
      bool alt=(i+j)%2; int t1[3]={a,b,alt?d:c}, t2[3]={alt?a:b,d,c};
      for(int*t: {t1,t2}){Tri tr; for(int k=0;k<3;++k)tr.n[k]=t[k]; tr.body=body; T.push_back(tr);} } };
  block(0,0,0); block(0.0+0.13,L+1e-9,1);   // B au-dessus de A, decale (contact oblique sur l'arete)
  int nn=X.size(); double Dm[3][3]; double f=E/((1+nu)*(1-2*nu)); 
  Eigen::Matrix3d D; D<<f*(1-nu),f*nu,0, f*nu,f*(1-nu),0, 0,0,f*(1-2*nu)/2;
  for(auto&t:T){ V2 P[3]={X[t.n[0]],X[t.n[1]],X[t.n[2]]}; double A2=pot::cross2(P[1]-P[0],P[2]-P[0]); t.A0=0.5*A2;
    double b[3],c[3]; for(int k=0;k<3;++k){int j=(k+1)%3,l=(k+2)%3; b[k]=P[j].y()-P[l].y(); c[k]=P[l].x()-P[j].x();}
    t.B.setZero(); for(int k=0;k<3;++k){t.B(0,2*k)=b[k]/A2; t.B(1,2*k+1)=c[k]/A2; t.B(2,2*k)=c[k]/A2; t.B(2,2*k+1)=b[k]/A2;}
    for(int k=0;k<3;++k) m[t.n[k]]+=rho*t.A0/3; }
  std::vector<V2> x=X, v(nn,V2(0,0)), F(nn);
  // chargement
  double sig=5e6; std::vector<V2> fext(nn,V2(0,0));
  int nA=(N+1)*(N+1);
  if(mode==0){ for(int i=nA;i<nn;++i) v[i]=V2(0,-2.0); }      // rebond : B tombe sur A a 2 m/s
  else { // presse : traction de compression sig sur la face basse de A et haute de B (forces nodales constantes)
    for(int i=0;i<=N;++i){ double w=(i==0||i==N)?0.5:1.0; fext[i]+=V2(0, sig*L/N*w); fext[nA+N*(N+1)+i]+=V2(0,-sig*L/N*w);} 
    std::mt19937 r2(3); std::normal_distribution<double> G(0,0.05); for(int i=0;i<nn;++i) v[i]=V2(G(r2),G(r2)); }
  // pas stable comme computeStableDt : 2 sqrt(m/(2E + 2 p)) * dtFactor, et CFL
  double dt=1e30; for(int i=0;i<nn;++i) dt=std::min(dt,2*std::sqrt(m[i]/(2*E+2*p))); 
  double cp=std::sqrt(E*(1-nu)/(rho*(1+nu)*(1-2*nu))); dt=std::min(dt,0.3*L/N/cp); dt*=dtfac;
  auto energies=[&](double&ke,double&ue,double&uc,double&wpot){ ke=ue=uc=0; for(int i=0;i<nn;++i) ke+=0.5*m[i]*v[i].squaredNorm();
    for(auto&t:T){ Eigen::Matrix<double,6,1> u; for(int k=0;k<3;++k){V2 d=x[t.n[k]]-X[t.n[k]]; u(2*k)=d.x();u(2*k+1)=d.y();} Eigen::Vector3d e=t.B*u; ue+=0.5*t.A0*e.dot(D*e);} 
    for(auto&a:T) if(a.body==0) for(auto&b:T) if(b.body==1){ V2 pa[3]={x[a.n[0]],x[a.n[1]],x[a.n[2]]},pb[3]={x[b.n[0]],x[b.n[1]],x[b.n[2]]}; double IA,IB,ar; if(enq::integrals(pa,pb,IA,IB,ar)) uc+=p*(IA+IB);} 
    wpot=0; for(int i=0;i<nn;++i) wpot-= fext[i].dot(x[i]-X[i]); };
  double ke,ue,uc,wp; energies(ke,ue,uc,wp); double E0=ke+ue+uc+wp; double Wc=0, Wx=0, bias=0;
  double Ttot = mode==0? 3e-3 : 0.05; long ns=(long)(Ttot/dt); double maxPen=0;
  for(long s=0;s<ns;++s){
    for(int i=0;i<nn;++i) F[i]=fext[i];
    for(auto&t:T){ Eigen::Matrix<double,6,1> u; for(int k=0;k<3;++k){V2 d=x[t.n[k]]-X[t.n[k]]; u(2*k)=d.x();u(2*k+1)=d.y();}
      Eigen::Matrix<double,6,1> fe=-t.A0*t.B.transpose()*(D*(t.B*u)); for(int k=0;k<3;++k) F[t.n[k]]+=V2(fe(2*k),fe(2*k+1)); }
    for(auto&a:T) if(a.body==0) for(auto&b:T) if(b.body==1){ V2 pa[3]={x[a.n[0]],x[a.n[1]],x[a.n[2]]},pb[3]={x[b.n[0]],x[b.n[1]],x[b.n[2]]};
      pot::PairForce R; if(!pot::pairForce(pa,pb,p,R)) continue; pot::PairForce R0=R; if(corr) enq::correct(pa,pb,p,R);
      for(int k=0;k<3;++k){ F[a.n[k]]+=R.fA[k]; F[b.n[k]]+=R.fB[k]; Wc+=(R.fA[k].dot(v[a.n[k]])+R.fB[k].dot(v[b.n[k]]))*dt;
        Wx+=((R0.fA[k]-R.fA[k]).dot(v[a.n[k]])+(R0.fB[k]-R.fB[k]).dot(v[b.n[k]]))*dt; } maxPen=std::max(maxPen,R.area); }
    for(int i=0;i<nn;++i){ bias+=F[i].squaredNorm()*dt*dt/(2*m[i]); v[i]+=dt/m[i]*F[i]; x[i]+=dt*v[i]; }
  }
  double ke1,ue1,uc1,wp1; energies(ke1,ue1,uc1,wp1); double E1=ke1+ue1+uc1+wp1;
  printf("mode=%d corr=%d p=%.1e dt=%.3e steps=%ld | E0=%.6e E1=%.6e dE=%+.4e (%.3e rel) | Wc=%+.4e Uc_fin=%.3e aire_max=%.2e Wnongrad=%+.3e\n",
     mode,corr,p,dt,ns,E0,E1,E1-E0,(E1-E0)/std::max(1e-30,std::fabs(mode==0?E0:ke+ue+uc)),Wc,uc1,maxPen,Wx);
}
