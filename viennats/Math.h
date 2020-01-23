#ifndef DEF_MYMATH_2D
#define DEF_MYMATH_2D

/* =========================================================================
   Copyright (c)    2008-2015, Institute for Microelectronics, TU Wien.

                            -----------------
                 ViennaTS - The Vienna Topography Simulator
                            -----------------

   Contact:         viennats@iue.tuwien.ac.at

   License:         MIT (X11), see file LICENSE in the base directory
============================================================================= */
#define DEBG 0

#include <cmath>
#include <array>
#include <chrono>
#include <random>
#include <fstream>
//at for fourRateInterpolation()
#include "LSlib/levelset.hpp"
//#include <boost/static_assert.hpp>

//delaunator
#include <delaunator.hpp>
#include <cstdio>
///Namespace for all custom mathematical and statistical tools.
namespace my {
  ///Includes mathematical functions like equation solvers, root finders and vector arithmetic.
  namespace math {

    const double Pi1_2   =1*std::asin(1.);
    const double Pi     =2*std::asin(1.);
    const double Pi2   =4*std::asin(1.);
    const double PiH2_8  =std::asin(1.)*std::asin(1.)*0.5;

    template<class T> int Sign(T);
    template<class T> int SignPos(T);
    template<class T> int SignNeg(T);
    template<class T> int SignEps(T,T);
    template<class T> T pow(T,int);
    template<class T> T pow2(T);
    int fac(int);
    int SignPow(int);

    int SolveCubic(double c[4], double s[3]);
    int SolveQudratic(double c[3],double s[2]);
    int SolveLinear(double c[2],double s[1]);

    ///factorial
    template<int n> class factorial {
    public:
      enum {RET = n * factorial<n - 1>::RET};
    };
    template<> class factorial<0>{
    public:
      enum { RET = 1 };
    };

    //template<int D>  class LinearRayInterpolation;
    template<class CoefT, int Order> class Polynom;

    double FindFirstRoot(double IntervalStart, double IntervalEnd, const Polynom<double, 2>& p);
    double FindFirstRoot(double IntervalStart, double IntervalEnd, const Polynom<double, 3>& p);

    template<int D> void TransformRho(double*);
    //template <int D> int GetRho(double* rho, double* distances, unsigned int* pts );

    template <int D> void DetermineCoefficientsForImplicitRayTracing(const double* Position, const double* Direction, const double* Rho, double* Coefficients );
    template <int D> void CalculateNormal(double* n, const double* v, const double* A);

    template<class T>  bool isEqual(T x, T y, T eps);

    /*template <int D, template <class, int> class V > inline V<int,D> round(const V<double,D>& vd) {
            V<int,D> vi;
            for (int i=0;i<D;i++) vi[i]=static_cast<int>(nearbyint(vd[i]));
            return vi;
        }

    template <int D, template <class, int> class V > inline V<int,D> floor(const V<double,D>& vd) {
      V<int,D> vi;
      for (int i=0;i<D;i++) vi[i]=static_cast<int>(std::floor(vd[i]));
      return vi;
    }

    template <int D, template <class, int> class V > inline V<int,D> ceil(const V<double,D>& vd) {
      V<int,D> vi;
      for (int i=0;i<D;i++) vi[i]=static_cast<int>(std::ceil(vd[i]));
      return vi;
    }*/

    template<class T>
    inline T Interpolate(const T* X, const T* Y, int N, T v)
    {

        if (v<=X[0]) return Y[0];
            if (v>=X[N-1]) return Y[N-1];

            int a=0;
            int b=N-1;

            while (a+1!=b) {
                int c=(a+b)/2;
                if (X[c]>v) b=c; else a=c;
            };

            return Y[a]+(Y[b]-Y[a])*((v-X[a])/(X[b]-X[a]));
        }



    template <int D>
    struct TransformRho2 {
        static void exec(double Rho[]) {
            enum {tmp=1<<(D-1)};
            TransformRho2<D-1>::exec(Rho);
            TransformRho2<D-1>::exec(Rho+tmp);
            for (int i=0;i<tmp;++i) Rho[i+tmp]-=Rho[i];
        }
    };

    template <>
        struct TransformRho2<0> {
            static void exec(double Rho[]) {}
        };
  }
}

template<class T> inline int my::math::Sign(T x) {
  if (x>T(0)) return 1;
  if (x<T(0)) return -1;
  return 0;
}

template<class T> inline int my::math::SignPos(T x) {
  if (x<T(0)) return -1;
  return 1;
}

template<class T> inline int my::math::SignNeg(T x) {
  if (x>T(0)) return 1;
  return -1;
}

template<class T> inline int my::math::SignEps(T x, T eps) {
  if (x>eps) {
    return 1;
  } else if (x<-eps) {
    return -1;
  } else {
    return 0;
  }
}


template<class T> inline T my::math::pow(T x, int k) {
  T p(1);
  for (int i=0;i<k;i++) p*=x;
  return p;
}

template<class T> inline T my::math::pow2(T x) {
  return x*x;
}


inline int my::math::fac(int k) {
  int p=1;
  for (int i=2;i<=k;i++) p*=i;
  return p;
}

inline int my::math::SignPow(int k) {
  if (k%2==0) return 1; else return -1;
}



inline int my::math::SolveCubic(double c[4], double s[3]) {
  if (std::fabs(c[3])<=1e-10*std::fabs(c[0])) {
    return SolveQudratic(c,s);
  } else {
    const double Pi2_3=(M_PI*2.)/3.;


    double A = c[ 2 ] / (3.*c[ 3 ]);
    double B = c[ 1 ] / c[ 3 ];
    double C = c[ 0 ] / c[ 3 ];

    double A2=A*A;
    double R=A2*A+(C-A*B)*0.5;
    double R2=R*R;
    double Q=(A2-B/3.);
    double Q3=Q*Q*Q;

    if (R2<Q3) {

      double Q1_2=sqrt(Q);
      double theta=acos(R/(Q*Q1_2))/3.;

      s[0]=-2.*Q1_2*cos(theta)-A;
      s[1]=-2.*Q1_2*cos(theta+Pi2_3)-A;
      s[2]=-2.*Q1_2*cos(theta-Pi2_3)-A;

      return 3;
    } else {

      double AA=-cbrt(R+my::math::Sign(R)*sqrt(R2-Q3));

      if (AA!=0.) {
        s[0]=AA+Q/AA-A;
        return 1;
      } else {
        s[0]=AA-A;
        return 1;
      }
    }
  }
}

int my::math::SolveQudratic(double c[3],double s[2]) {
  if (std::fabs(c[2])<=1e-10*std::fabs(c[0])) {
    return SolveLinear(c,s);
  } else {
    double D=c[1]*c[1]-4.*c[0]*c[2];
    if (D<0.) {
      return 0;
    } else {
      double Q=-0.5*(c[1]+my::math::Sign(c[1])*std::sqrt(D));
      s[0]=Q/c[2];
      s[1]=c[0]/Q;

      return 2;
    }
  }
}

int my::math::SolveLinear(double c[2],double s[1]) {
  if (std::fabs(c[1])<=1e-10*std::fabs(c[0])) {
    return 0;
  } else {
    s[0]=-c[0]/c[1];
    return 1;
  }
}

template<class T> inline
bool  my::math::isEqual(T x, T y, T eps){
    return std::fabs(x-y) <= eps;
}



namespace my {
  namespace math {

    /*template <int D> LinearRayInterpoaltion {

    public:

      LinearRayInterpolation()



    }*/
    template<class CoefT, int Order>  class Polynom {
      CoefT coefficients[Order+1];
    public:
      const CoefT* Coefficients() const {
        return coefficients;
      }

      CoefT* Coefficients() {
        return coefficients;
      }

      CoefT Coefficients(int order) const {
          return coefficients[order];
      }

      Polynom() {}

      template <class ArgT> CoefT operator()(ArgT x) const {
        CoefT result=coefficients[Order];
        for (int g=Order-1;g>=0;g--) {
          result*=x;
          result+=coefficients[g];
        }
        return result;
      }

      template <class ArgT> CoefT Derivate1(ArgT x) const {
        CoefT result=Order*coefficients[Order];
        for (int g=Order-1;g>0;g--) {
          result*=x;
          result+=g*coefficients[g];
        }
        return result;
      }
    };


    template <class T> T FindFirstTransitionFromPosToNegOfPolynomBisect(T t0, T t1, const Polynom<T, 3>& p ,T eps) {

      if (t0>t1) return std::numeric_limits<T>::max();

      T f=p(t0);

      if (f<-eps) {
        return t0;
      } else if (f<=eps) {
        if (p.Derivate1(t0)<=0.) return t0;
      }

      f=p(t1);

      bool Intersection=false;
      if (f<-eps) {
        Intersection=true;
      } else if (f<=eps) {
        if (p.Derivate1(t1)<=0.) {
          Intersection=true;
        }
      }

      T Discriminant=p.Coefficients()[2]*p.Coefficients()[2]-3.*p.Coefficients()[1]*p.Coefficients()[3];
      if (Discriminant>=0) {  //p' has real roots

        //calculate roots of p'
        T Q=-(p.Coefficients()[2]+my::math::SignPos(p.Coefficients()[2])*std::sqrt(Discriminant));

        T e0=Q/(3.*p.Coefficients()[3]);

        T e1=p.Coefficients()[1]/Q;

        if (e0>e1) std::swap(e0,e1);

        if ((e0<t1) && (e0>t0)) { //if e0 element of [t0,t1]
          if (p(e0)>0.) {
            t0=e0;
          } else {
            t1=e0;
            Intersection=true;
          }
        }

        if ((e1<t1) && (e1>t0)) { //if e1 element of [t0,t1]
          if (p(e1)>0.) {
            t0=e1;
          } else {
            t1=e1;
            Intersection=true;
          }
        }

      }

      if (!Intersection) return std::numeric_limits<T>::max();

      while(std::fabs(t1-t0)>eps) {
        T t=0.5*(t1+t0);
        f=p(t);
        if (f==0.) {
          return t;
        } else if (f>0.) {
          t0=t;
        } else {
          t1=t;
        }
      }

      return t0;

    }



    template <class T> T FindFirstTransitionFromPosToNegOfPolynomNewton(T t0, T t1, const Polynom<T, 3>& p ,T eps) {

      if (t0>t1) return std::numeric_limits<T>::max();

      T f=p(t0);

      if (f<-eps) {
        return t0;
      } else if (f<=eps) {
        if (p.Derivate1(t0)<=0.) return t0;
      }

      f=p(t1);

      bool Intersection=false;
      if (f<-eps) {
        Intersection=true;
      } else if (f<=eps) {
        if (p.Derivate1(t1)<=0.) {
          Intersection=true;
        }
      }

      T Discriminant=p.Coefficients(2)*p.Coefficients(2)-3.*p.Coefficients(1)*p.Coefficients(3);
      if (Discriminant>=0) {  //p' has real roots

        //calculate roots of p'
        T Q=-(p.Coefficients(2)+my::math::SignPos(p.Coefficients(2))*std::sqrt(Discriminant));

        T e0=Q/(3.*p.Coefficients(3));

        T e1=p.Coefficients()[1]/Q;

        if (e0>e1) std::swap(e0,e1);

        if ((e0<t1) && (e0>t0)) { //if e0 element of [t0,t1]
          if (p(e0)>0.) {
            t0=e0;
          } else {
            t1=e0;
            Intersection=true;
          }
        }

        if ((e1<t1) && (e1>t0)) { //if e1 element of [t0,t1]
          if (p(e1)>0.) {
            t0=e1;
          } else {
            t1=e1;
            Intersection=true;
          }
        }

      }

      if (!Intersection) return std::numeric_limits<T>::max();

      T t=(t1+t0)*0.5;
      T dxold=t1-t0;
      T dx=dxold;

      T ft=p(t);
      T dft=p.Derivate1(t);

      while(true) {

        if ((   ((t-t0)*dft-ft)*((t-t1)*dft-ft)>0.   ) || (std::fabs(ft+ft)>std::fabs(dxold*dft))) {
          dxold=dx;
          dx=0.5*(t1-t0);
          t=t0+dx;
          if ((t0==t) || (dx<eps)) return t0;
        } else {
          dxold=dx;
          dx=ft/dft;
          T temp=t;
          t-=dx;
          if ((temp==t) || (std::fabs(dx)<eps)) return std::max(t-eps,t0);
        }


        ft=p(t);
        dft=p.Derivate1(t);

        if (ft>0.0) {
          t0=t;
        } else {
          t1=t;
        }
      }

    }

    template <class T> T FindFirstTransitionFromPosToNegOfPolynomNewton(T t0, T t1, const Polynom<T, 2>& p ,T eps) {

      if (t0>t1) return std::numeric_limits<T>::max();

            T f=p(t0);

            if (f<-eps) {
                return t0;
            } else if (f<=eps) {
                if (p.Derivate1(t0)<=0.) return t0;
            }

            f=p(t1);

            bool Intersection=false;
            if (f<-eps) {
                Intersection=true;
            } else if (f<=eps) {
                if (p.Derivate1(t1)<=0.) {
                    Intersection=true;
                }
            }

            T e=-0.5*(p.Coefficients(1)/p.Coefficients(2));        //determine inflection point

            if ((e<t1) && (e>t0)) { //if e0 element of [t0,t1]
                if (p(e)>0.) {
                    t0=e;
                } else {
                    t1=e;
                    Intersection=true;
                }
            }

            if (!Intersection) return std::numeric_limits<T>::max();

            T t=(t1+t0)*0.5;
            T dxold=t1-t0;
            T dx=dxold;

            T ft=p(t);
            T dft=p.Derivate1(t);

            while(true) {

                if ((   ((t-t0)*dft-ft)*((t-t1)*dft-ft)>0.   ) || (std::fabs(ft+ft)>std::fabs(dxold*dft))) {
                    dxold=dx;
                    dx=0.5*(t1-t0);
                    t=t0+dx;
                    if ((t0==t) || (dx<eps)) return t0;
                } else {
                    dxold=dx;
                    dx=ft/dft;
                    T temp=t;
                    t-=dx;
                    if ((temp==t) || (std::fabs(dx)<eps)) return std::max(t-eps,t0);
                }


                ft=p(t);
                dft=p.Derivate1(t);

                if (ft>0.0) {
                    t0=t;
                } else {
                    t1=t;
                }
            }
    }

    double FindFirstRoot(double t0, double t1, const Polynom<double, 3>& p) {

      //BOOST_STATIC_ASSERT(std::numeric_limits<double>::has_infinity);

      const int N=5;

      double f0=p(t0);
      int sgn0=Sign(f0);

      if (sgn0<0) return t0;

      double f1=p(t1);
      int sgn1=Sign(f1);

      double Discriminant=p.Coefficients()[2]*p.Coefficients()[2]-3.*p.Coefficients()[1]*p.Coefficients()[3];
      if (Discriminant>=0) {  //p' has real roots

        //calculate roots of p'
        double Q=-(p.Coefficients()[2]+my::math::SignPos(p.Coefficients()[2])*std::sqrt(Discriminant));

        double e0=Q/(3.*p.Coefficients()[3]);

        double e1=p.Coefficients()[1]/Q;

        if (e0>e1) std::swap(e0,e1);

        if ((e0<=t1) && (e0>=t0)) { //if e0 element of [t0,t1]
          double fe0=p(e0);
          int sgn_fe0=Sign(fe0);
          if (sgn_fe0==sgn0) {
            t0=e0;
            f0=fe0;
            sgn0=sgn_fe0;
          } else {
            t1=e0;
            f1=fe0;
            sgn1=sgn_fe0;
          }
        }

        if ((e1<=t1) && (e1>=t0)) { //if e1 element of [t0,t1]
          double fe1=p(e1);
          int sgn_fe1=Sign(fe1);
          if (sgn_fe1==sgn0) {
            t0=e1;
            f0=fe1;
            sgn0=sgn_fe1;
          } else {
            t1=e1;
            f1=fe1;
            sgn1=sgn_fe1;
          }
        }

      }


      if (sgn0==0) {
        return t0;
      } else if (sgn1==0) {
        return t1;
      } else if (sgn0==sgn1) {
        return std::numeric_limits<double>::max();
      }

      for (int i=0;i<N;i++) {
        //double t=t0-(t1-t0)*(f0/(f1-f0));
        double t=(t1+t0)/2;
        double ft=p(t);
        int sgn_ft=Sign(ft);

        if (sgn_ft==0) return t;

        if (sgn_ft==sgn0) {
          t0=t;
          f0=ft;
        } else {
          t1=t;
          f1=ft;
        }
      }
      return t0-(t1-t0)*(f0/(f1-f0));
      //return (t1+t0)/2;



    }

    double FindFirstRoot(double t0, double t1, const Polynom<double, 2>& p) {

      //BOOST_STATIC_ASSERT(std::numeric_limits<double>::has_infinity);

      double Discriminant=p.Coefficients()[1]*p.Coefficients()[1]-4.*p.Coefficients()[0]*p.Coefficients()[2];
      if (Discriminant>=0) {  //p has real roots

        //calculate roots of p
        double Q=-(p.Coefficients()[1]+my::math::SignPos(p.Coefficients()[1])*std::sqrt(Discriminant))*0.5;

        double e0=Q/p.Coefficients()[2];

        double e1=p.Coefficients()[0]/Q;

        if (e0>e1) std::swap(e0,e1);

        if ((e0<=t1) && (e0>=t0)) { //if e0 element of [t0,t1]
          return e0;
        }

        if ((e1<=t1) && (e1>=t0)) { //if e1 element of [t0,t1]
          return e1;
        }


      }
      return std::numeric_limits<double>::max();
    }


    template <> void DetermineCoefficientsForImplicitRayTracing<2>(const double* Position, const double* Direction, const double* Rho, double* Coefficients ) {

      Coefficients[2]=Direction[0]*Rho[3];
      Coefficients[1]=Position[1]*Coefficients[2];
      Coefficients[2]*=Direction[1];
      Coefficients[0]=Position[0]*Rho[3]+Rho[2];
      Coefficients[1]+=Direction[1]*Coefficients[0]+Direction[0]*Rho[1];
            Coefficients[0]*=Position[1];
            Coefficients[0]+=Position[0]*Rho[1]+Rho[0];

            /*
            //TEST
            double t=43.3;
            double x=Position[0]+t*Direction[0];
            double y=Position[1]+t*Direction[1];
            double val1=Coefficients[2]*t*t+Coefficients[1]*t+Coefficients[0];
            double val2=x*y*Rho[3]+y*Rho[2]+x*Rho[1]+Rho[0];
            assert(std::fabs(val1-val2)<1e-8);
            */

    }








    /*template <> void TransformRho<3>(double* Rho) {

      double g=Rho[0];
      g-=Rho[1];
      g-=Rho[2];
      g-=Rho[4];

      Rho[7]-=g;
      Rho[7]-=Rho[3];
      Rho[7]-=Rho[5];
      Rho[7]-=Rho[6];

      Rho[3]+=g;
      Rho[3]+=Rho[4];

      Rho[5]+=g;
      Rho[5]+=Rho[2];

      Rho[6]+=g;
      Rho[6]+=Rho[1];

      Rho[1]-=Rho[0];
      Rho[2]-=Rho[0];
      Rho[4]-=Rho[0];

    }

    template <> void TransformRho<2>(double* Rho) {
      Rho[1]-=Rho[0];
      Rho[3]-=Rho[1];
      Rho[3]-=Rho[2];
      Rho[2]-=Rho[0];
    }*/

    template <> void DetermineCoefficientsForImplicitRayTracing<3>(const double* Position, const double* Direction, const double* Rho, double* Coefficients ) {

      Coefficients[2]=Direction[1]*Direction[0];
      Coefficients[3]=Coefficients[2]*Direction[2]*Rho[7];
      Coefficients[2]*=(Position[2]*Rho[7]+Rho[3]);

      double tmp=Position[1]*Rho[7]+Rho[5];

      Coefficients[0]=tmp*Position[0]+Position[1]*Rho[6]+Rho[4];
      Coefficients[1]=Coefficients[0]*Direction[2];
      Coefficients[0]*=Position[2];

      tmp*=Direction[0];
      tmp+=(Direction[1]*(Position[0]*Rho[7]+Rho[6]));

      Coefficients[2]+=Direction[2]*tmp;
      Coefficients[1]+=Position[2]*tmp;

      tmp=Position[0]*Rho[3]+Rho[2];

      Coefficients[1]+=((Position[1]*Rho[3]+Rho[1])*Direction[0]+Direction[1]*tmp);
      Coefficients[0]+=(Position[1]*tmp+Rho[1]*Position[0]+Rho[0]);

    }

    template <> void CalculateNormal<3>(double* n, const double* v, const double* Rho) {
      n[0]=v[1]*Rho[7]+Rho[5];
      n[2]=n[0];

      n[0]*=v[2];
      n[0]+=(v[1]*Rho[3]+Rho[1]);

      n[1]=v[2]*(v[0]*Rho[7]+Rho[6])+v[0]*Rho[3]+Rho[2];

      n[2]*=v[0];
      n[2]+=(v[1]*Rho[6]+Rho[4]);

      double no2=n[0]*n[0]+n[1]*n[1]+n[2]*n[2];
      if (no2>0.) {
        no2=std::sqrt(no2);
        n[0]/=no2;
        n[1]/=no2;
        n[2]/=no2;
      } else {
        n[0]=0.;
        n[1]=0.;
        n[2]=0.;
      }
    }

    template <> void CalculateNormal<2>(double* n, const double* v, const double* Rho) {
      n[0]=Rho[1]+v[1]*Rho[3];
      n[1]=Rho[2]+v[0]*Rho[3];

      double no2=n[0]*n[0]+n[1]*n[1];

      if (no2>0.) {
                no2=std::sqrt(no2);
                n[0]/=no2;
                n[1]/=no2;
            } else {
                n[0]=0.;
                n[1]=0.;
            }
    }

    //Four rate Hubbart Interpolation. Algorithm works only in 3D, input 2D vector (vx,vy) is treated as (vx,vy,0)
    template<class T,int D>
    T fourRateInterpolation(lvlset::vec<T,D> nv, lvlset::vec<T,3> direction100, lvlset::vec<T,3> direction010, T r100, T r110, T r111, T r311 ){

        T Velocity=0;
        lvlset::vec<T,3> directions[3];

        directions[0] = direction100;
        directions[1] = direction010;

        directions[0]=Normalize(directions[0]);
        directions[1]=Normalize(directions[1]-directions[0]*dot(directions[0],directions[1]));
        directions[2]=cross(directions[0], directions[1]);

        lvlset::vec<T,3> NormalVector;
        NormalVector[0] = nv[0];
        NormalVector[1] = nv[1];
        if(D==3){
          NormalVector[2] = nv[2];
        }else{
          NormalVector[2] = 0;
        }

        NormalVector=Normalize(NormalVector);

        for (int i=0;i<3;++i) assert(dot(directions[i], directions[(i+1)%3])<1e-6);

        lvlset::vec<T,3> N;

        for (int i=0;i<3;i++) N[i]=std::fabs(directions[i][0]*NormalVector[0]+directions[i][1]*NormalVector[1]+directions[i][2]*NormalVector[2]);
        N.reverse_sort();

        assert(std::fabs(Norm(N)-1)<1e-4);


        if (dot(N, lvlset::vec<T,3>(-1,1,2))<0) {
            Velocity=-((r100*(N[0]-N[1]-2*N[2])+r110*(N[1]-N[2])+3*r311*N[2])/N[0]);    //region A
        } else {
            Velocity=-((r111*((N[1]-N[0])*0.5+N[2])+r110*(N[1]-N[2])+1.5*r311*(N[0]-N[1]))/N[0]);//region C
        }

        return Velocity;
    }

    //Weighted essentially non-oscillatory differentiation scheme 3rd order
    //x1 ... x5 stencil points from left to right
    //plus == true => right-sided
    template<class T>
    inline T weno3(T x1, T x2, T x3, T x4, T x5, T dx, bool plus=true, T eps=1e-6){

      T dxp1 = x2 - x1; //i-2
      T dxp2 = x3 - x2; //i-1
      T dxp3 = x4 - x3; //i
      T dxp4 = x5 - x4;  //i+1


      T result = 0;

      if(plus==true){

          T rp = (eps + math::pow2( dxp4 - dxp3 )) / (eps + math::pow2(dxp3 - dxp2));
          T wp = 1.0 / (1 + 2.0 * math::pow2(rp));

          result = dxp2 + dxp3 - wp * (dxp4 - 2.0 * dxp3 + dxp2);

          result /= (2.0 * dx);

      } else {
          T rp = (eps + math::pow2( dxp2 - dxp1)) / (eps + math::pow2( dxp3 - dxp2 ));
          T wp = 1.0 / (1 + 2.0 * math::pow2(rp));

          result = dxp2 + dxp3 - wp * (dxp1 - 2.0 * dxp2 + dxp3);

          result /= (2.0 * dx);
      }

      return result;
    }

    //Weighted essentially non-oscillatory differentiation scheme 5th order
    //x1 ... x7 stencil points from left to right
    //plus == true => right-sided
    template<class T>
    inline T weno5(T x1, T x2, T x3, T x4, T x5, T x6, T x7, T dx, bool plus=true, T eps=1e-6){

      if(plus==false){
        T v1 = (x2 - x1)/dx; //i-3
        T v2 = (x3 - x2)/dx; //i-2
        T v3 = (x4 - x3)/dx; //i-1
        T v4 = (x5 - x4)/dx; //i
        T v5 = (x6 - x5)/dx;  //i+1


        T p1 = v1/3.0 - 7*v2/6.0 + 11*v3/6.0;
        T p2 = -v2/6.0 + 5*v3/6.0 + v4/3.0;
        T p3 = v3/3.0 + 5*v4/6.0 - v5/6.0;

        T s1 = 13/12.0  * pow2(v1 - 2*v2 + v3) + 1/4.0 * pow2(v1 - 4*v2 + 3*v3);
        T s2 = 13/12.0  * pow2(v2 - 2*v3 + v4) + 1/4.0 * pow2(v2 - v4);
        T s3 = 13/12.0  * pow2(v3 - 2*v4 + v5) + 1/4.0 * pow2(3*v3 - 4*v4 + v5);

        T al1 = 0.1/(eps + s1);
        T al2 = 0.6/(eps + s2);
        T al3 = 0.3/(eps + s3);

        T alsum = al1 + al2 + al3;

        T w1=al1 / alsum;
        T w2=al2 / alsum;
        T w3=al3 / alsum;

        return w1 * p1 + w2 * p2 + w3 * p3;
      } else{
        T v1 = (x7 - x6)/dx;
        T v2 = (x6 - x5)/dx;
        T v3 = (x5 - x4)/dx;
        T v4 = (x4 - x3)/dx;
        T v5 = (x3 - x2)/dx;


        T p1 = v1/3.0 - 7*v2/6.0 + 11*v3/6.0;
        T p2 = -v2/6.0 + 5*v3/6.0 + v4/3.0;
        T p3 = v3/3.0 + 5*v4/6.0 - v5/6.0;

        T s1 = 13/12.0  * pow2(v1 - 2*v2 + v3) + 1/4.0 * pow2(v1 - 4*v2 + 3*v3);
        T s2 = 13/12.0  * pow2(v2 - 2*v3 + v4) + 1/4.0 * pow2(v2 - v4);
        T s3 = 13/12.0  * pow2(v3 - 2*v4 + v5) + 1/4.0 * pow2(3*v3 - 4*v4 + v5);

        T al1 = 0.1/(eps + s1);
        T al2 = 0.6/(eps + s2);
        T al3 = 0.3/(eps + s3);

        T alsum = al1 + al2 + al3;

        T w1=al1 / alsum;
        T w2=al2 / alsum;
        T w3=al3 / alsum;

        return w1 * p1 + w2 * p2 + w3 * p3;

      }
    }



    //returns spherical vector [r,theta,phi], with theta: elevation, phi: azimuth
    template<class T> inline
    lvlset::vec<T,3> cartesianToSpherical(const lvlset::vec<T,3>& cart){
      T r = Norm(cart);

      T theta = acos(cart[2]/r);
      T phi = atan2(cart[1],cart[0]);
      if(phi < 0)
        phi += 2*math::Pi;

      return lvlset::vec<T,3>{r,theta,phi};
    }

    //takes spherical vector [r,theta,phi], with theta: elevation, phi: azimuth
    template<class T> inline
    lvlset::vec<T,3> sphericalToCartesian(const lvlset::vec<T,3>& sph){

      return lvlset::vec<T,3>{    sph[0]*sin(sph[1])*cos(sph[2]),
                                  sph[0]*sin(sph[1])*sin(sph[2]),
                                  sph[0]*cos(sph[1])};
    }



    //check if a point is on spherical triangle that is located on a unit sphere.
    //The triangle is given by its vertices (cartesian coordinates). Orientation of vertices has to be right-handed!
    template<class T>
    bool isOnSphericalTriangle(const lvlset::vec<T,3>& point,
                               const lvlset::vec<T,3>& vertex1,
                               const lvlset::vec<T,3>& vertex2,
                               const lvlset::vec<T,3>& vertex3){

          const T eps=1e-4;
          T distances[3] = {0, 0, 0};
          const lvlset::vec<T,3> vertices[3] = {vertex1,vertex2, vertex3};

          //distance = v cdot (pi x p(i+1)) for the triangle vertices pi and point v.
          //If point v is on the same side of all 3 planes defined by (0, pi, p(i+1)),
          //then it is on the spherical triangle.

          //Distance loop
          for(int j = 0 ; j < 3; ++j){
               //Cross product loop
               for(int i = 0 ; i < 3; ++i){
                   distances[j] += vertices[(0+j)%3][(1+i)%3]*vertices[(1+j)%3][(2+i)%3]*point[(0+i)%3] -
                                   vertices[(0+j)%3][(2+i)%3]*vertices[(1+j)%3][(1+i)%3]*point[(0+i)%3];
               }
          }

         // if(isEqual(point(0),point(),1e-12))
       //    std::cout << "isOnSphericalTriangle: point:" << point << ", distances: " << distances[0] << ", " << distances[1] << ", " << distances[2] << std::endl;


          if ( (distances[0] >= 0 || isEqual(distances[0],0.0,eps)) &&
               (distances[1] >= 0 || isEqual(distances[1],0.0,eps)) &&
               (distances[2] >= 0 || isEqual(distances[2],0.0,eps))){
              return true;
          }

          return false;
     }


    template<class T>
    lvlset::vec<T,3> sphericalBarycentricCoords(const lvlset::vec<T,3>& v, const lvlset::vec<T,3>& v1, const lvlset::vec<T,3>& v2, const lvlset::vec<T,3>& v3 ){

       lvlset::vec<T,3> normals[3]{cross(v2,v3), cross(v3,v1), cross(v1,v2)};
       lvlset::vec<T,3> sin_al;

       for(size_t i = 0; i < 3; ++i)
          sin_al[i] = dot(v,normals[i]);

       lvlset::vec<T,3> sin_be{dot(v1,normals[0]), dot(v2,normals[1]), dot(v3,normals[2]) };

       lvlset::vec<T,3> result;
       for(size_t i = 0; i < 3; ++i)
          result[i] = sin_al[i] / sin_be[i];

       return result;
     }

    //Stereographic projection, projection pole is SOUTH pole
    template<class T>
    lvlset::vec<T,2> stereographicToPlane(const lvlset::vec<T,3>& v){
        return lvlset::vec<T,2>{2*v[0] / (1 + v[2]),2*v[1] / (1 + v[2]) };
    }

    template<class T>
    lvlset::vec<T,3> stereographicToSphere(const lvlset::vec<T,2>& v){
        T s = 4/( v[0]*v[0] + v[1]*v[1] + 4);
        return lvlset::vec<T,3>{s*v[0], s*v[1], 2*s - 1};
    }

  }

  namespace symmetry {


    //Convert Miller-Bravais indices (4 entries) to Cartesian vector
    //a1 ... base vector in c plane
    //c ... base vector along c direction
    //NOTE: Magnitude of a1 and c define c/a ratio of hexagonal/trigonal crystal system
    template<class T>
    lvlset::vec<T,3> millerBravaisToCartesian(const std::array<T,4>& hex, const lvlset::vec<T,3>& a1, const lvlset::vec<T,3>& c){

      lvlset::vec<T,3> a2, a3;

      if(std::fabs(hex[0] + hex[1] + hex[2]) > 1e-4){
        std::cerr << "millerBravaisToCartesian: invalid input [" << hex[0] << ", " << hex[1] << ", " << hex[2] << ", " << hex[3] << "]\n";
        exit(-1);
      }

      a2 = lvlset::RotateAroundAxis(a1, c, 2.0 * math::Pi/3);
      a3 = lvlset::RotateAroundAxis(a1, c, 4.0 * math::Pi/3);

      return hex[0] * a1 + hex[1] * a2 + hex[2] * a3 + hex[3] * c;
    }

    //calculate normal on a plane (Miller Bravais) given by hex for a hexagonal system with lattice parameter ratio c/a = ca_ratio
    template<class T>
    std::array<T,4> millerBravaisNormalVector(const std::array<T,4>& hex, const T ca_ratio){
      T lastEntry = 1.5*hex[3]/(ca_ratio*ca_ratio);

     return std::array<T,4>{hex[0],hex[1],-hex[0]-hex[1],lastEntry}; 
    }


    //Trigonal symmetry for Sapphire
    template <class T> class D3d {

      private:
         const lvlset::vec<T,3> c3_1{0,0,1}; //3 fold rotation
         const lvlset::vec<T,3> a1{sqrt(3)*0.5, -0.5, 0};
         const lvlset::vec<T,3> sigma1 = a1;


         //NB. angle between a1 and sigma1 is 90 deg (ensured by constructor)

         const T c3_1_angle = 2*math::Pi/3;

         //INT_NUM + 1, due to [1,0,-1,0] and [1,-1,0,0] being equivalent.
         //Both directions are required to describe the entire fundmental domain,
         // but interpolation values have to be equal
         std::vector<std::array<T,4>> interpPlanesHex;
         const std::array<T,4> MM_plane{1,0,-1,0};//additional m direction 
          
//          { {
//                                                                  {0,0,0,1},//c direction
//                                                                  {1,-1,0,2},//r direction
//                                                                  {1,-1,0,0},//m direction
//                                                                  {1,-1,0,5},//Shen 1 direction
//                                                                  {4,-5,1,38},//Shen 3 direction
//                                                                  {1,-1,0,12},//Shen 4 direction
//                                                                  {1,0,-1,5},//11 direction
//                                                                  {1,-1,0,38},//O direction
//                                                                  {1,1,-2,38},//V direction
//                                                                  {1,0,-1,0}//additional m direction 
//                                                                    } };
//

         std::vector<std::array<size_t,3>> interpSphTri;
       //  { {
       //                                                            {1, 2, 9},
       //                                                            {2, 6, 9},
       //                                                            {3, 6, 2},
       //                                                            {3, 4, 6},
       //                                                            {5, 4, 3},
       //                                                            {6, 7, 0},
       //                                                            {4, 7, 6},
       //                                                            {5, 7, 4},
       //                                                            {7, 8, 0},
       //                                                            {8, 7, 5},
       //                                                            } };
         std::vector<std::array<T,4>> interpDirectionsHex;
         std::vector<lvlset::vec<T,3>> interpVertices;
         std::vector<std::vector<T>> interpRates; //rate vector for every material
         
         T basalAngle = 0; //angle between user defined a1 and a1{sqrt(3)*0.5, -0.5, 0}
         T ca_ratio = 1; //ratio between lattice parameters c and a, this depends on the material


         //sampling arrays
         //plain velocity
         using SamplingTable = std::vector<std::vector<T> >; 

         std::vector<SamplingTable> vel_sampled;
         std::vector<SamplingTable> vel100_sampled; //n + dN*ex
         std::vector<SamplingTable> velm100_sampled; //n - dN*ex
         std::vector<SamplingTable> vel010_sampled; //n + dN*ey
         std::vector<SamplingTable> vel0m10_sampled; //n -dN*ey
         std::vector<SamplingTable> vel001_sampled; //n +dN*ez
         std::vector<SamplingTable> vel00m1_sampled; //n -dN*ez

         const T phi_min = 0;
         const T phi_max = 2*math::Pi;

         const T u_min = -1;
         const T u_max = 1;

         size_t sampleM = 0;//default value
         T du = 0;//default value
         T dphi = 0;//default value
        
         const T CEPS = std::cbrt(std::numeric_limits<T>::epsilon());
      public:

          //in current status, defineCoordinateSystem HAS TO BE called in order to guanrantee meaningful state of object
          D3d(){
          }

          //a ... normalized vector in c plane
          //ca ... c/a lattice parameter ratio 
          void defineCoordinateSystem(const lvlset::vec<T,3>&  a, const T ca ){

              basalAngle = lvlset::SignedAngle(a, a1,c3_1);
              ca_ratio = ca;
              std::cout << "a (user specified system)= " << a << ", a1 (internal system)= " << a1 << ", c3_1 (internal system)= " << c3_1 << ", Basal angle = " << basalAngle << ", c/a = " << ca_ratio << "\n";

              
            }

          void defineRateFunction(const std::vector<std::vector<T>>& planes, const std::vector<std::vector<T>>& rates){
             
             size_t mplane_index=planes.size()+1; 

              for(size_t i = 0; i < planes.size(); ++i){
                  if(planes[i] == std::vector<T>{1, -1, 0, 0})
                      mplane_index=i;
                      
                  interpPlanesHex.push_back(std::array<T,4>{planes[i][0],planes[i][1],planes[i][2],planes[i][3]});
              }

              if(mplane_index > planes.size()){
                  std::cout << "Error: M plane not given\n";
                  exit(0);
              }



              interpPlanesHex.push_back(MM_plane);//M' plane is required

              //TODO check for C and M plane

              interpRates = rates;
              for(size_t m=0; m < rates.size(); ++m){
                  interpRates[m].push_back(rates[m][mplane_index]);
              }


              //Set interpolation vectors again with new coordinate system 
              for(size_t i=0; i < interpPlanesHex.size() ; ++i){
                  interpDirectionsHex.push_back( millerBravaisNormalVector(interpPlanesHex[i], ca_ratio));
                  lvlset::vec<T,3> vec =  lvlset::Normalize(millerBravaisToCartesian(interpDirectionsHex[i],a1,ca_ratio * c3_1));
                  interpVertices.push_back(reduceToFundmental(vec));
              }

              // std::cout << "Angle([1 -1 0 5], [0 0 0 1]) = Angle(" << interpVertices[4] << ", " << interpVertices[0] << ") = " << Angle(interpVertices[4], interpVertices[0]) << "\n";

              //Triangulate interpVertices
              //
              //Stereographic projection
              /* x0, y0, x1, y1, ... */
              std::vector<T> stereoPlaneVertices;
              
              for(size_t i=0; i < interpVertices.size(); ++i)
              {
                  auto sp =  my::math::stereographicToPlane(interpVertices[i]);
                  stereoPlaneVertices.push_back(sp[0]);
                  stereoPlaneVertices.push_back(sp[1]);
              }

              std::cout << "stereoPlaneVertices:\n";
              for( auto s : stereoPlaneVertices){
                  std::cout << s << "\n";
              }
              
               //triangulation happens here
              delaunator::Delaunator d(stereoPlaneVertices);
              for(std::size_t i = 0; i < d.triangles.size(); i+=3) {
                          printf(
                              "Triangle points: [[%f, %f], [%f, %f], [%f, %f]]\n",
                              d.coords[2 * d.triangles[i]],        //tx0
                              d.coords[2 * d.triangles[i] + 1],    //ty0
                              d.coords[2 * d.triangles[i + 1]],    //tx1
                              d.coords[2 * d.triangles[i + 1] + 1],//ty1
                              d.coords[2 * d.triangles[i + 2]],    //tx2
                              d.coords[2 * d.triangles[i + 2] + 1] //ty2
                              );

                          
                              }

              //build interpSphTri based on Delaunator result
              for(size_t i = 0; i < d.triangles.size(); i+=3){
                  T tx0 =  d.coords[2 * d.triangles[i]];        //tx0
                  T ty0 = d.coords[2 * d.triangles[i] + 1];    //ty0
                  T tx1 = d.coords[2 * d.triangles[i + 1]];    //tx1
                  T ty1 =d.coords[2 * d.triangles[i + 1] + 1];//ty1
                  T tx2 = d.coords[2 * d.triangles[i + 2]];    //tx2
                  T ty2 = d.coords[2 * d.triangles[i + 2] + 1];//ty2

                  T area = 0.5*std::fabs( tx0*ty1 - tx0*ty2 - tx1*ty0 + tx1*ty2 + tx2*ty0 -tx2*ty1 ); 
                  std::cout << "area = " << area << "\n";

                  if(area > 1e-6){
                      interpSphTri.push_back(std::array<size_t,3>{d.triangles[i],d.triangles[i+1], d.triangles[i+2]});
                  }
              }

              std::cout << "size interpSphTri = " << interpSphTri.size() << ", size interpVertices = " << interpVertices.size() << "\n";

              std::cout << "interpSphTri:\n";

              for(auto in : interpSphTri){
                  for(auto e : in){
                      std::cout << e << " ";
                  }
                  std::cout << "\n";
              }
              std::cout << "\n";

              std::cout << "interpRates[0]:\n";
              for(auto r : interpRates[0])
                  std::cout << r << " ";

              std::cout << "\n";


          }

          void sampleRateFunctions(const size_t M) {


                size_t mats(interpRates.size());

                du = (u_max-u_min)/M;
                dphi = (phi_max-phi_min)/M;

                auto t1=std::chrono::system_clock::now();

                std::cout << "Started sampling of D3d rate function, M = " << M << "\n";
                for(size_t matNum=0; matNum < mats; ++matNum){

                    //TODO more elegant sampling table 
                    SamplingTable vel(M+1);
                    for(auto& v : vel) v.resize(M+1);

                    SamplingTable vel100(M+1);
                    for(auto& v : vel100) v.resize(M+1);

                    SamplingTable velm100(M+1);
                    for(auto& v : velm100) v.resize(M+1);
                    
                    SamplingTable vel010(M+1);
                    for(auto& v : vel010) v.resize(M+1);

                    SamplingTable vel0m10(M+1);
                    for(auto& v : vel0m10) v.resize(M+1);

                    SamplingTable vel001(M+1);
                    for(auto& v : vel001) v.resize(M+1);

                    SamplingTable vel00m1(M+1);
                    for(auto& v : vel00m1) v.resize(M+1);


                    for(size_t i=0; i <= M; ++i){                   
                        T phi = phi_min + i*dphi;
                        
                        for(size_t j=0; j <= M; ++j){
                            T u = u_min + j*du;

                            T nx = sqrt(1-u*u)*cos(phi);
                            T ny = sqrt(1-u*u)*sin(phi);
                            T nz = u; 

                            T v =  interpolate(lvlset::vec<T,3>{nx,ny,nz},matNum);                         
                            T DN = v * CEPS;

                            vel[i][j] = v;
                           
                             
                            vel100[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx+DN,ny,nz}), matNum);
                                        

                            
                            velm100[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx-DN,ny,nz}), matNum);
    

                            vel010[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx,ny+DN,nz}), matNum);
                            
                            vel0m10[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx,ny-DN,nz}), matNum);
                            
                            vel001[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx,ny,nz+DN}), matNum);
                            
                            vel00m1[i][j] =  interpolate(Normalize(lvlset::vec<T,3>{nx,ny,nz-DN}), matNum);
                        }

                    }

                  vel_sampled.push_back(vel); 
                  vel100_sampled.push_back(vel100);
                  velm100_sampled.push_back(velm100);
                  vel010_sampled.push_back(vel010);
                  vel0m10_sampled.push_back(vel0m10);
                  vel001_sampled.push_back(vel001);
                  vel00m1_sampled.push_back(vel00m1);

                  std::cout << "Material #" << matNum << " done\n";
                      

                }
                sampleM = M;
                std::cout << "Finished D3d sampling\n";

                auto t2=std::chrono::system_clock::now();
                std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
                std::cout << "Sampling took " << fp_ms.count() << " ms\n";
           


                int N = 100000;
                std::default_random_engine generator;
                generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
                std::uniform_real_distribution<T> dist_phi(0,2*math::Pi);
                std::uniform_real_distribution<T> dist_u(-1,1);

                std::vector<lvlset::vec<T,3>> nvec(N);
                std::vector<T> interpolationResult(N);
 
                for(int i(0); i < N; ++i){
                    T phi=dist_phi(generator);
                    T u=dist_u(generator);
 
                    T nx = sqrt(1-u*u) * cos(phi);
                    T ny = sqrt(1-u*u) * sin(phi);
                    T nz = u;
 
                    nvec[i] = lvlset::vec<T,3>{nx,ny,nz};
                }
                for(int j(0); j < N; ++j){
                    interpolationResult[j] = interpolateSampled(nvec[j], 0);
                }
                

                ofstream file("rateFunction.XYZ");
                file << "nx,ny,nz,interpolationResult\n";
                if(file.is_open()){
                    for(int j(0); j < N; ++j){
                        file << nvec[j][0] << "," << nvec[j][1] << "," << nvec[j][2] << "," << interpolationResult[j]<< "\n";
                    }
 
                }
           }
 
        //test function to time the effect of sampling
           void timingTestSampling(const int N, 
                                  const int matNum){
               
               std::vector<lvlset::vec<T,3>> nvec(N);
               std::vector<T> interpolationResult(N);
               std::vector<T> samplingResult(N);

               std::cout << "Sampling test\n";

               std::default_random_engine generator;
               generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
               std::uniform_real_distribution<T> dist_phi(0,2*math::Pi);
               std::uniform_real_distribution<T> dist_u(-1,1);

               for(int i(0); i < N; ++i){
                   T phi=dist_phi(generator);
                   T u=dist_u(generator);

                   T nx = sqrt(1-u*u) * cos(phi);
                   T ny = sqrt(1-u*u) * sin(phi);
                   T nz = u;

                   nvec[i] = lvlset::vec<T,3>{nx,ny,nz};
               }

               std::cout << "\tRNG finshed\n\tStarting default interpolation run\n"; 
               auto t1=std::chrono::system_clock::now();

               for(int j(0); j < N; ++j){
                   interpolationResult[j] = interpolate(nvec[j], matNum);
               }
               
               auto t2=std::chrono::system_clock::now();
               std::chrono::duration<double, std::milli> fp_ms = t2 - t1;
               std::cout << "\tDefault interpolation took " << fp_ms.count() << " ms\n";

               std::cout << "\tDefault interpolation run finshed\n\tStarting sampled access run\n"; 
               auto t3=std::chrono::system_clock::now();

               for(int j(0); j < N; ++j){
                   samplingResult[j] = interpolateSampled(nvec[j], matNum);
               }
               
               auto t4=std::chrono::system_clock::now();
               std::chrono::duration<double, std::milli> fp_ms2 = t4 - t3;
               std::cout << "\tSampled access took " << fp_ms2.count() << " ms\n";

               ofstream file("samplingTestResults.csv");
               file << "nx, ny, nz, interpolationResult, samplingResult, absDiff, relDiff\n";
               if(file.is_open()){
                   for(int j(0); j < N; ++j){
                       file << nvec[j][0] << ", " << nvec[j][1] << ", " << nvec[j][2] << ", " << interpolationResult[j] << ", " << samplingResult[j] << ", " << fabs(interpolationResult[j] - samplingResult[j]) << ", " << fabs( (interpolationResult[j] - samplingResult[j])/interpolationResult[j]) << "\n";
                   }

               }


           } 
          
          //reduce to fundamental domain, which is 0<theta<Pi/2, 0<phi<2*PI/3
           lvlset::vec<T,3> reduceToFundmental(lvlset::vec<T,3> in) const {
              lvlset::vec<T,3> out = in;
              lvlset::vec<T,3> out_sph = math::cartesianToSpherical(in);

              if(out[2] < 0){
                out = -in;
                out_sph = math::cartesianToSpherical(out);
              }

              while(out_sph[2] > 2.0*math::Pi/3){
                out = lvlset::RotateAroundAxis(out,c3_1, -c3_1_angle);
                out_sph = math::cartesianToSpherical(out);
              }

              if(out_sph[2] > math::Pi/3.0){
                out = lvlset::ReflectionInPlane(out,sigma1);
              }

              return out;
           }

           T interpolateSampled(const lvlset::vec<T,3> in, const int materialNum) const {
                T u = in[2];
                T phi = atan2(in[1],in[0]);
                if(phi < 0)
                   phi += 2*math::Pi;  

                int i = static_cast<int>( std::floor((phi-phi_min)/dphi + 0.5));
                int j = static_cast<int>( std::floor((u-u_min)/du + 0.5));

                return vel_sampled[materialNum][i][j];
           }

           T interpolateSLFSampled(const lvlset::vec<T,3> in, const int materialNum, const int ix, const int iy, const int iz) const {
                lvlset::vec<T,3> in_norm = Normalize(in);
                T u = in_norm[2];
                T phi = atan2(in_norm[1],in_norm[0]);
                if(phi < 0)
                   phi += 2*math::Pi;  

                int i = static_cast<int>( std::floor((phi-phi_min)/dphi + 0.5));
                int j = static_cast<int>( std::floor((u-u_min)/du + 0.5));

                if( ix == 1 && iy == 0 && iz == 0)
                    return vel100_sampled[materialNum][i][j];
                
                if( ix == -1 && iy == 0 && iz == 0)
                    return velm100_sampled[materialNum][i][j];

                if( ix == 0 && iy == 1 && iz == 0)
                    return vel010_sampled[materialNum][i][j];
                
                if( ix == 0 && iy == -1 && iz == 0)
                    return vel0m10_sampled[materialNum][i][j];
                
                if( ix == 0 && iy == 0 && iz == 1)
                    return vel001_sampled[materialNum][i][j];
                
                if( ix == 0 && iy == 0 && iz == -1)
                    return vel00m1_sampled[materialNum][i][j];

                std::cerr << "interpolateSLFSampled, invalid (ix, iy, iz) combination: " << ix << ", " << iy << ", " << iz << "\n";
                exit(-1);
                return 0;
           }


           //input vector is assumed to be normalized |v|=1
           T interpolate(const lvlset::vec<T,3> in, const int materialNum) const{

              lvlset::vec<T,3> in_fund  = lvlset::RotateAroundAxis(in, c3_1, -basalAngle);
              in_fund = reduceToFundmental(in_fund);

#if 0
              if(dot(in_fund,lvlset::vec<T,3>{0,0,1})/Norm(in_fund)  > 0.993){
                  return rC;
              }

#endif
              bool triangleFound=false;
              size_t triangleIdx = 0;

              for( ; triangleIdx < interpSphTri.size(); ++triangleIdx){
                if(math::isOnSphericalTriangle(in_fund, interpVertices[interpSphTri[triangleIdx][2]],
                                                        interpVertices[interpSphTri[triangleIdx][1]],
                                                        interpVertices[interpSphTri[triangleIdx][0]])){

                    triangleFound=true;
                    break;
                }
              }

              if(!triangleFound){
                std::cerr << "\nTriangle not found\n";
                std::cerr << "Point (fundmental): " << in_fund << "\n";
                exit(-1);
              }

              lvlset::vec<T,3> baryCoords = math::sphericalBarycentricCoords(in_fund, interpVertices[interpSphTri[triangleIdx][2]],
                                                                       interpVertices[interpSphTri[triangleIdx][1]],
                                                                       interpVertices[interpSphTri[triangleIdx][0]]);
              T result = 0;
              T sum = 0;

              for(size_t i=0; i<3; ++i){
                result+=baryCoords[i] * interpRates[materialNum][ interpSphTri[triangleIdx][2-i] ]; //linear interpolation
                sum += baryCoords[i];
              }
              result /= sum; //we divide by  b0 + b1 + b2  to ensure partition of unity property

              return result;
           }

    };
  }
}






#endif //DEF_MYMATH_2D
