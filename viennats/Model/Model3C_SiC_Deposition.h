    /*
     * ModelSapphireWetEtching.h
     *
     *  Created on: 1/2020
     *      Author: atoifl
     */

#ifndef MODEL3C_SIC_DEPOSITION_H_
#define MODEL3C_SIC_DEPOSITION_H_

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
            const bool use_sampling=true;
            const long M = 200; //samples in phi and cos(theta)

            //Sapphire lattice parameters
            const double a = 4.785;//Angstrom
            const double c = 12.991;//Angstrom
            
            lvlset::vec<double,3> directionA;
            lvlset::vec<double,3> directionC;

            std::vector<std::vector<double>> planes; //planes given by experiment
            std::vector<std::vector<double>> rates; //one vector per material number, stores rates for planes given by experiment
            
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
                
                std::cout << Parameters;

                bool b = parse(
                        Parameters.begin(),
                        Parameters.end(),
                        *(
                              (str_p("directionA")  >> '='  >> '{' >> real_p[assign_a(directionA[0])]  >> "," >> real_p[assign_a(directionA[1])] >> "," >> real_p[assign_a(directionA[2])] >> '}' >> ';') |
                              (str_p("directionC")  >> '='  >> '{' >> real_p[assign_a(directionC[0])]  >> "," >> real_p[assign_a(directionC[1])] >> "," >> real_p[assign_a(directionC[2])] >> '}' >> ';') 
                        ),
                        space_p | comment_p("//") | comment_p("/*", "*/")).full;

                //if (!b) msg::print_error("Failed interpreting process parameters!");

                //Parse planes= and rates= using a dirty handwritten parser
                std::vector<std::vector<double>> planes;
                const int EntryNum = 4;
                std::string planeStr="planes=";
                parseListOfList(Parameters,planeStr,EntryNum,planes);
               
                std::cout << "Parsed planes:\n"; 
                for(auto p : planes){
                    for(auto v : p){
                        std::cout << v << " ";
                    }

                    std::cout << "\n";
                }


                std::vector<std::vector<double>> rates;
                int ratesEntryNum = planes.size();
                std::string rateStr="rates=";
                parseListOfList(Parameters,rateStr,ratesEntryNum,rates);
                
                std::cout << "Parsed rates:\n"; 
                for(auto p : rates){
                    for(auto v : p){
                        std::cout << v << " ";
                    }

                    std::cout << "\n";
                }


                sapphireSymmetry.defineCoordinateSystem(directionA,c/a);
                sapphireSymmetry.defineRateFunction(planes,rates);
                
                if(use_sampling){
                  sapphireSymmetry.sampleRateFunctions(M);

                    if(false)
                        sapphireSymmetry.timingTestSampling(10000000, 0 );
                }

                //reverse because material id is reversed w.r.t. layer id (for output)
                /*std::reverse(r100.begin(),r100.end());
                std::reverse(r110.begin(),r110.end());
                std::reverse(r111.begin(),r111.end());
            std::reverse(r311.begin(),r311.end());*/

            // find materials with no growth in any direction and store in zeroVel
#if 0
            for(unsigned int i=0; i < rC.size(); ++i){
                zeroVel.push_back(false);
              if(fabs(rC[i]) < EPS)
                if(fabs(rR[i]) < EPS)
                  if(fabs(rM[i]) < EPS)
                      if(fabs(rO[i]) < EPS)
                          if(fabs(r1[i]) < EPS)
                            if(fabs(r3[i]) < EPS)
                              if(fabs(r4[i]) < EPS)
                                if(fabs(r11[i]) < EPS)
                                    if(fabs(rV[i]) < EPS)
                                        zeroVel[i]=true;
            }
#endif

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

#if 0
            if (zeroVel[Material] == true) {
                Velocity=0;
                return;
            }
#endif
            lvlset::vec<double,3> nv{NormalVector[0],NormalVector[1],NormalVector[2]};

            //Velocity using standard interpolation
            if(!use_sampling)
                Velocity = -sapphireSymmetry.interpolate(nv,Material);
            else
                Velocity = -sapphireSymmetry.interpolateSampled(nv,Material);
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

       template<class VecType>
        void CalculateSLFVelocity(
                double &Velocity,
                const VecType& NormalVector,
                const double *Coverages,
                const double *Rates,
                int Material,
                bool connected,
                bool visible,
                const int ix, const int iy, const int iz) const {

            if(!use_sampling){
                std::cerr << "Error: trying to use sampled rates, but sampling is turned off in ModelSapphire\n";
                exit(-1);
            }
#if 0
            if (zeroVel[Material] == true) {
                Velocity=0;
                return;
            }
#endif
            lvlset::vec<double,3> nv{NormalVector[0],NormalVector[1],NormalVector[2]};

            Velocity = -sapphireSymmetry.interpolateSLFSampled(nv,Material,ix,iy,iz);

        }


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

        void parseListOfList(const std::string in, const std::string quantityName, const int EntryNum, std::vector<std::vector<double>>& parsed){
        

           
            std::cout << " ---- '" << in << "' ----\n";
            auto pos1 = in.find(quantityName);
            auto pos2 = in.find(";", pos1);

            std::string input = in.substr(pos1+quantityName.size(),pos2-pos1-quantityName.size());
            
                
            std::cout << " ---- '" << input << "' ----\n";

            int state = 0;

            for(std::size_t i=0; i < input.size(); ++i){
                
                
                if(state == 2)
                { 
                   std::vector<double> arr(EntryNum); 
                   for(std::size_t j = 0; j < EntryNum; ++j){
                        
                       auto commaPos = (j != EntryNum - 1) ? input.find(",",i) : input.find("}",i);
                       auto sub = input.substr(i,commaPos-i);
                       std::cout << sub << "\n";
                       double entry = std::stod(sub);
                       arr[j]=entry;
                       i = commaPos + 1;
                   }

                   parsed.push_back(arr);
                   state = 1;               
                }

                
                if(state ==1 && input[i] == '{')
                    state = 2;
                if(state == 0 && input[i] == '{')
                   state = 1;
            }
        }

    };

//    const unsigned int SapphireWetEtching::NumberOfParticleClusters[ConstantRates::NumberOfParticleTypes]={};

}
#endif /* MODEL3C_SIC_DEPOSITION_H__ */
