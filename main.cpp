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
   Trotter::init(tau);

   PEPS<double> peps;
   peps.initialize_state(D);

   peps.normalize(D_aux);

   Environment::calc_env('A',peps,D_aux);

   cout << Heisenberg::energy(peps) << endl;

   DArray<2> Sz(d,d);
   Sz(0,0) = -0.5;
   Sz(0,1) = 0.0;
   Sz(1,0) = 0.0;
   Sz(1,1) = 0.5;

   for(int i = 0;i < 500;++i){

      propagate::step(peps,D_aux);

      peps.normalize(D_aux);

      Environment::calc_env('A',peps,D_aux);

      cout << i << "\t" << Heisenberg::energy(peps) << "\t" << Heisenberg::local(peps,Sz) << endl;

   }

   Environment::calc_env('A',peps,D_aux);
   Environment::test_env();

}
