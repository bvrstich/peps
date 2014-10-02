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

   //initialize some statics dimensions
   Global::lat.set(L,L,d);
   Environment::init();
   Heisenberg::init();

   for(int f_i = 7000;f_i < 8000;++f_i){

      double f = f_i / 10000.0;

      PEPS<double> peps;

      peps.set_jastrow(f);
      peps.normalize(D_aux);

      Environment::calc_env('A',peps,D_aux);

      cout << f << "\t" << Heisenberg::energy(peps)/(double)(L*L) << endl;

   }

}
