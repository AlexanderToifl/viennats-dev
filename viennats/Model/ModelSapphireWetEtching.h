/*
 * ModelSapphireWetEtching.h
 *
 *  Created on: 9/2019
 *      Author: atoifl
 */

#ifndef MODELSAPPHIREWETETCHING_H_
#define MODELSAPPHIREWETETCHING_H_

#include <stack>
#include <vector>
#include <cassert>
//#include <algorithm>
//#include <functional>

#include <boost/spirit/include/classic.hpp>
#include "../message.h"
#include "../parser_actors.h"
#include "vector.hpp"

namespace model {

///Anisotropic wet etching model for sapphire etching using H2SO4 : H3PO4

    class SapphireWetEtching {

        const double EPS = 1e-6;

        lvlset::vec<double,3> directionA;
        lvlset::vec<double,3> directionC;

        //[1,0,-1,0], [0,0,0,1], [1,-1,0,5], [4,-5,1,38], [1,-1,0,12], [1,-1,0,0]
        std::vector<double> r10m10;
        std::vector<double> r0001;
        std::vector<double> r1m105;
        std::vector<double> r4m5138;
        std::vector<double> r1m1012;

        std::vector<bool> zeroVel;

        my::symmetry::D3d<double> sapphireSymmetry;


    public:

        static const bool OutputFluxes=false;
        static const bool SpatiallyEqualDistributedFlux=true;
        static const bool ReemissionIsMaterialDependent=true;	//TODO is set to true for output
        static const bool CalculateConnectivities=false;        //TODO can be useful?
        static const bool CalculateVisibilities=false;
        static const bool CalculateNormalVectors=true;

        class ParticleType {
        public:
            double Direction[3];
            double Flux;
        };

        SapphireWetEtching(const std::string & Parameters) {
            using namespace boost::spirit::classic;
            using namespace parser_actors;


            bool b = parse(
                    Parameters.begin(),
                    Parameters.end(),
                    *(
                          (str_p("directionA")  >> '='  >> '{' >> real_p[assign_a(directionA[0])]  >> "," >> real_p[assign_a(directionA[1])] >> "," >> real_p[assign_a(directionA[2])] >> '}' >> ';') |
                          (str_p("directionC")  >> '='  >> '{' >> real_p[assign_a(directionC[0])]  >> "," >> real_p[assign_a(directionC[1])] >> "," >> real_p[assign_a(directionC[2])] >> '}' >> ';') |
                          (str_p("rate10m10")  >>  '='  >>  '{' >> ( real_p[push_back_a(r10m10)]  % ',')>> '}'  >> ';') |
                          (str_p("rate0001")  >>  '='  >>  '{' >> ( real_p[push_back_a(r0001)]  % ',')>> '}'  >> ';') |
                          (str_p("rate1m105")  >>  '='  >>  '{' >> ( real_p[push_back_a(r1m105)]  % ',')>> '}'  >> ';') |
                          (str_p("rate4m5138")  >>  '='  >>  '{' >> ( real_p[push_back_a(r4m5138)]  % ',')>> '}'  >> ';') |
                          (str_p("rate1m1012")  >>  '='  >>  '{' >> ( real_p[push_back_a(r1m1012)]  % ',')>> '}'  >> ';')
                    ),
                    space_p | comment_p("//") | comment_p("/*", "*/")).full;

            if (!b) msg::print_error("Failed interpreting process parameters!");

            sapphireSymmetry.initInterpolationTriangles(directionA, directionC);

            //reverse because material id is reversed w.r.t. layer id (for output)
            /*std::reverse(r100.begin(),r100.end());
            std::reverse(r110.begin(),r110.end());
            std::reverse(r111.begin(),r111.end());
            std::reverse(r311.begin(),r311.end());*/

            // find materials with no growth in any direction and store in zeroVel
            for(unsigned int i=0; i < r10m10.size(); ++i){
              zeroVel.push_back(false);
              if(fabs(r10m10[i]) < EPS)
                if(fabs(r0001[i]) < EPS)
                  if(fabs(r1m105[i]) < EPS)
                    if(fabs(r4m5138[i]) < EPS)
                      if(fabs(r1m1012[i]) < EPS)
                        zeroVel[i]=true;
            }


        }

        static const int CoverageStorageSize=0;
        static const int RatesStorageSize=0;
        static const unsigned int NumberOfParticleTypes=0;
        unsigned int NumberOfParticleClusters[1];

        template<class VecType>
        void CalculateVelocity(
                double &Velocity,
                const VecType& NormalVector,
                const double *Coverages,
                const double *Rates,
                int Material,
                bool connected,
                bool visible) const {


            if (zeroVel[Material] == true) {
                Velocity=0;
                return;
            }

            lvlset::vec<double,3> nv{NormalVector[0],NormalVector[1],NormalVector[2]};

            // lvlset::vec<T,3> NormalVector;
            // NormalVector[0] = nv[0];
            // NormalVector[1] = nv[1];
            // if(D==3){
            //   NormalVector[2] = nv[2];
            // }else{
            //   NormalVector[2] = 0;
            // }

            Velocity = -sapphireSymmetry.interpolate(nv,r10m10[Material], r0001[Material], r1m105[Material], r4m5138[Material],  r1m1012[Material]);

            //Velocity = my::symmetry::sapphireFiveRateInterpolation<double,3>(nv, directionA, directionC, r10m10[Material], r0001[Material], r1m105[Material], r4m5138[Material],  r1m1012[Material]);

            //std::cout << Velocity << "\n";

        }

		template<class VecType>
		void CalculateVectorVelocity(
				double *Velocity,
				const VecType& NormalVector,
				const double *Coverages,
				const double *Rates,
				int Material,
				bool connected,
				bool visible) const {}


        static void UpdateCoverage(double *Coverages, const double *Rates) {}

        template <class PT> static void ParticleGeneration(PT& p, int ParticleType, double ProcessTime, double* Position) {}

        template <class PT, class NormVecType> void ParticleCollision(
                                    const PT& p,
                                    const NormVecType& NormalVector,
                                    double* Rates,
                                    const double* Coverages,
                                    double RelTime) const {}

		template <class PT, class VecType> static void ParticleReflexion(
							const PT& p,
							std::stack<PT>& particle_stack,
							const VecType& NormalVector,
							const double* Coverages,
							int Material//,
//                            int D,
//                            double dot // dot product between the incoming particle direction and the normal vector
                            ) {}
    };

//    const unsigned int SapphireWetEtching::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}
#endif /* MODELSAPPHIREWETETCHING_H_ */
