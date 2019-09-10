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


#include <cmath>
#include <array>


//at for fourRateInterpolation()
#include "LSlib/levelset.hpp"
//#include <boost/static_assert.hpp>

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
    //The triangle is given by its vertices (cartesian coordinates).
    template<class T>
    bool isOnSphericalTriangle(const lvlset::vec<T,3>& point,
                               const lvlset::vec<T,3>& vertex1,
                               const lvlset::vec<T,3>& vertex2,
                               const lvlset::vec<T,3>& vertex3){

          const T eps=1e-12;
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

         // if(isEqual(point(0),point(1),1e-12))
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


  }

  namespace symmetry {

    //c-axis orientation and a1 direction orientation can be given, a2 and a3 are rotated by 120 deg around in c-axis
    //Default values lead to following setup:
    //x-axis: [1 0 -1 0]
    //y-axis: [-1 2 -1 0]
    //z-axis: [0 0 0 1]
    template<class T>
    //lvlset::vec<T,3> millerBravaisToCartesian(const std::array<T,4>& hex, lvlset::vec<T,3>& a1 = {sqrt(3)*0.5, -0.5, 0}, lvlset::vec<T,3>& c = {0, 0, 1} ){
    lvlset::vec<T,3> millerBravaisToCartesian(const std::array<T,4>& hex, const lvlset::vec<T,3>& a1, const lvlset::vec<T,3>& c){

      lvlset::vec<T,3> a2, a3;

      assert(std::fabs(hex[0] + hex[1] + hex[2]) > 1e-4);

      a2 = lvlset::RotateAroundAxis(a1, c, 2.0 * math::Pi/3);
      a3 = lvlset::RotateAroundAxis(a1, c, 4.0 * math::Pi/3);

      return hex[0] * a1 + hex[1] * a2 + hex[2] * a3 + hex[3] * c;
    }

    //Trigonal symmetry for Sapphire
    template <class T> class D3d {

      private:
         static constexpr size_t INT_NUM = 5;
         const lvlset::vec<T,3> c3_1; //3 fold rotation
         const lvlset::vec<T,3> a1;
         lvlset::vec<T,3> sigma1;

         //NB. angle between a1 and sigma1 is 90 deg (ensured by constructor)

         const T c3_1_angle = 2*math::Pi/3;

         //INT_NUM + 1, due to [1,0,-1,0] and [1,-1,0,0] being equivalent.
         //Both directions are required to describe the entire fundmental domain,
         // but interpolation values have to be equal
         const std::array<std::array<T,4>, INT_NUM+1> interpDirectionsHex{ {{1,0,-1,0}, //=r10m10
                                                                  {0,0,0, 1}, //=r0001
                                                                  {1,-1,0,5}, //=r1m105
                                                                  {4,-5,1,38},//=r4m5138
                                                                  {1,-1,0,12},//=r1m1012
                                                                  {1,-1,0,0}} };//=== r10m10

         const std::array<std::array<size_t,3>, INT_NUM> interpSphTri{ {{0,5,2},
                                                   {0,2,3},
                                                   {3,2,4},
                                                   {3,4,1},
                                                   {0,3,1}} };

         lvlset::vec<T,3> interpVertices[INT_NUM+1];


        public:
        //constructor
        //Default values result in
        // x-axis: [1 0 -1 0]
        // y-axis: [-1 2 -1 0]
        // z-axis: [0 0 0 1]
        //and mirror plane along [0 1 -1 0]
        D3d(  const lvlset::vec<T,3>&  a = {sqrt(3)*0.5, -0.5, 0} , const lvlset::vec<T,3>&  c = {0,0,1}) : c3_1(c), a1(a){

          sigma1 = RotateAroundAxis(a,c,math::Pi/2);



          for(size_t i=0; i < INT_NUM+1; ++i){

            interpVertices[i] = lvlset::Normalize(millerBravaisToCartesian(interpDirectionsHex[i],a1,c3_1));
            interpVertices[i] = reduceToFundmental(interpVertices[i]);
          }

        }

          //reduce to fundamental domain, which is 0<theta<Pi/2, 0<phi<2*PI/3
           lvlset::vec<T,3> reduceToFundmental(lvlset::vec<T,3> in){
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

           //input vector is assumed to be normalized |v|=1
           T interpolate(const lvlset::vec<T,3> in, T r10m10, T r0001, T r1m105, T r4m5138, T r1m1012){

              T interpValues[INT_NUM + 1] = {r10m10, r0001, r1m105, r4m5138, r1m1012,r10m10};
              lvlset::vec<T,3> in_fund = reduceToFundmental(in);

              bool triangleFound=false;
              size_t triangleIdx = 0;

              for( ; triangleIdx < INT_NUM + 1; ++triangleIdx){

                if(math::isOnSphericalTriangle(in_fund, interpVertices[interpSphTri[triangleIdx][0]],
                                                        interpVertices[interpSphTri[triangleIdx][1]],
                                                        interpVertices[interpSphTri[triangleIdx][2]]))
                  triangleFound=true;
                  break;
              }

              assert(triangleFound);

              lvlset::vec<T,3> baryCoords = math::sphericalBarycentricCoords(in_fund, interpVertices[interpSphTri[triangleIdx][0]],
                                                                       interpVertices[interpSphTri[triangleIdx][1]],
                                                                       interpVertices[interpSphTri[triangleIdx][2]]);

              T result = 0;
              T sum = 0;

              for(size_t i=0; i<3; ++i){
                result+=baryCoords[i] * interpValues[ interpSphTri[triangleIdx][i] ]; //linear interpolation
                sum += baryCoords[i];
              }
              result /= sum; //we divide by  b0 + b1 + b2  to ensure partition of unity property

              return result;
           }

    };

    //Five rate interpolation for saphhire etching
    template<class T,int D>
    T sapphireFiveRateInterpolation(const lvlset::vec<T,D>& nv, const lvlset::vec<T,3>& directionA, const lvlset::vec<T,3>& directionC,
                                    T r10m10, T r0001, T r1m105, T r4m5138, T r1m1012 ){

      D3d<T> sapphire{directionA, directionC};

      lvlset::vec<T,3> NormalVector;
      NormalVector[0] = nv[0];
      NormalVector[1] = nv[1];
      if(D==3){
        NormalVector[2] = nv[2];
      }else{
        NormalVector[2] = 0;
      }

      return sapphire.interpolate(NormalVector,r10m10,r0001,r1m105,r4m5138,r1m1012);

    }
  }
}






#endif //DEF_MYMATH_2D
