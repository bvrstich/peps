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

   double tau = 0.001;

   //initialize some statics dimensions
   global::init(D,D_aux,d,L,L,tau);

   double f = 0.74;

   PEPS<double> peps;
   peps.initialize_jastrow(f);
   peps.normalize();

   propagate::step(false,peps);

   return 0;
}
