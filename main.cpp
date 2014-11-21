/**
 * @mainpage 
 * This is an implementation of an imaginary time evolution algorithm for the optimization of PEPS.
 * @author Brecht Verstichel
 * @date 25-03-2014
 */
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <complex>

using std::cout;
using std::endl;
using std::vector;
using std::complex;
using std::ofstream;

#include "include.h"

using namespace btas;

int main(int argc,char *argv[]){

   cout.precision(15);

   int L = atoi(argv[1]);//dimension of the lattice: LxL
   int d = atoi(argv[2]);//physical dimension
   int D = atoi(argv[3]);//virtual dimension
   int D_aux = atoi(argv[4]);//virtual dimension

   bool update = true;
   double tau = 0.01;
   int n_steps = 10;

   int mult = D_aux / (D*D);

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   PEPS<double> peps;
   peps.initialize_ising(0,0.001);

   for(int i= 0;i < 2000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 2000;i < 10000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 10000;i < 20000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   //increase D
   D++;
   peps.grow_bond_dimension(D,0.001);

   D_aux = mult * D * D;

   tau = 0.01;

   global::init(D,D_aux,d,L,L,tau);

   for(int i= 20000;i < 22000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 22000;i < 30000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }

   tau /= 10.0;
   global::stau(tau);

   for(int i= 30000;i < 40000;++i){

      propagate::step(update,peps,10);

      if(i % 10 == 0){

         peps.normalize();

         global::env.calc('A',peps);
         cout << i << "\t" << peps.energy()/(L*L) << endl;

      }

   }


   return 0;

}
