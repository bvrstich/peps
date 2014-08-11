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
   int D_aux = atoi(argv[4]);//auxiliary dimension for the contraction

   double tau = 0.01;

   //initialize some statics dimensions
   Global::lat.set(L,L,d);
   Environment::init();
   Heisenberg::init();

   Trotter::set(tau);

   PEPS<double> peps(D);
   peps.initialize_state_sum(D);

   peps.normalize(D_aux);

   char filename[200];

   sprintf(filename,"output/%dx%d/D=%d",L,L,D);

   //for(int i = 0;i < 1000;++i){

      //propagate::step(peps,D_aux);

      //if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << 0 << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      //}

   //}

/*
   tau /= 10.0;

   Trotter::set(tau);

   cout << endl;
   cout << "now for dt = " << tau << endl;
   cout << endl;

   for(int i = 1000;i < 2000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.01;
   Trotter::set(tau);
   peps.grow_bond_dimension(++D,0.001);

   sprintf(filename,"output/%dx%d/D=%d",L,L,D);

   D_aux = D*D;

   for(int i = 2000;i < 3000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.001;
   Trotter::set(tau);

   for(int i = 3000;i < 4000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.01;
   Trotter::set(tau);
   peps.grow_bond_dimension(++D,0.001);

   D_aux = D*D;

   for(int i = 4000;i < 5000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.001;
   Trotter::set(tau);

   for(int i = 5000;i < 6000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.01;
   Trotter::set(tau);
   peps.grow_bond_dimension(++D,0.001);

   D_aux = D*D;

   for(int i = 6000;i < 7000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }

   tau = 0.001;
   Trotter::set(tau);

   for(int i = 7000;i < 8000;++i){

      propagate::step(peps,D_aux);

      if(i % 10 == 0){

         Environment::calc_env('A',peps,D_aux);
         cout << i << "\t" << Heisenberg::energy(peps) << endl;

         //save:
         peps.save(filename);

      }

   }
*/
}
